# load libraries
library(forecast)
library(tidyverse)
library(tseries)
library(ggplot2)

# load inflation data
cpi_target_full <- read_csv("../data/cpi_target_full.csv")

# convert to time series data
y_ts <- ts(cpi_target_full$CPI_t, start = c(1961, 1), frequency = 12)
y_ts <- na.omit(y_ts)
summary(y_ts)

# plot inflation rate across time
autoplot(y_ts) +
  labs(title = "Transformed Inflation Rate (CPI_t)", 
       x = "Year", y = "Inflation Rate") +
  theme_minimal()

# check stationarity and variance
adf.test(y_ts)

# calculate rmse
rmse <- function(pred, truth) {
  sqrt(mean((truth - pred)^2))
}

# --------- AR(P) MODEL WITH FIXED TRAIN/TEST SPLIT ---------
# create lagged features
max_lag <- 12
lagged <- embed(y_ts, max_lag + 1)
targets <- lagged[, 1]
lag_matrix <- lagged[, -1, drop = FALSE]
colnames(lag_matrix) <- paste0("lag", seq_len(max_lag))

# split train/test
# last 66 months as test set, specifically from Jan 2020 onwards
nprev <- 100
total_rows <- nrow(lagged)
test_rows <- seq(total_rows - nprev + 1, total_rows)
train_rows <- seq_len(total_rows - nprev)

# double check date alignment
tail(cpi_target_full$sasdate, nprev)

# storing of results
results <- data.frame(
  p = integer(),
  AIC = numeric(),
  BIC = numeric(),
  HQ = numeric(),
  RMSE = numeric()
)

for (p in seq_len(max_lag)) {
  lag_subset <- lag_matrix[, seq_len(p), drop = FALSE]
  train_df <- data.frame(
    y = targets[train_rows],
    lag_subset[train_rows, , drop = FALSE]
  )
  model <- lm(y ~ ., data = train_df)

  test_df <- data.frame(
    y = targets[test_rows],
    lag_subset[test_rows, , drop = FALSE]
  )
  preds <- predict(model, newdata = test_df)

  rmse_val <- rmse(preds, test_df$y)
  aic_val <- AIC(model)
  bic_val <- BIC(model)

  n_train <- nrow(train_df)
  rss <- sum(resid(model)^2)
  sigma2 <- rss / n_train
  k_param <- length(coef(model))
  hq_val <- log(sigma2) + (2 * k_param * log(log(n_train))) / n_train

  results <- rbind(
    results,
    data.frame(
      p = p,
      AIC = aic_val,
      BIC = bic_val,
      HQ = hq_val,
      RMSE = rmse_val
    )
  )
}

best_aic <- results$p[which.min(results$AIC)]
best_bic <- results$p[which.min(results$BIC)]
best_hq <- results$p[which.min(results$HQ)]
best_rmse <- results$p[which.min(results$RMSE)]

cat(
  "Observations used:", length(y_ts), "\n",
  "Forecast horizon (nprev):", nprev, "\n\n"
)

print(results)

cat(
  "\nBest lag by AIC:", best_aic,
  "\nBest lag by BIC:", best_bic,
  "\nBest lag by HQ:", best_hq,
  "\nBest lag by out-of-sample RMSE:", best_rmse,
  "\n"
)

# plot lag vs ic
# reshape to long format for easier plotting
results_long <- tidyr::pivot_longer(
  results,
  cols = c("AIC", "BIC", "HQ", "RMSE"),
  names_to = "Criterion",
  values_to = "Value"
)

# plot AIC/BIC/HQ
ggplot(results_long %>% filter(Criterion %in% c("AIC", "BIC")),
       aes(x = p, y = Value, color = Criterion)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Model Selection Criteria by Lag Length (AR(p))",
    x = "Lag Order (p)",
    y = "Criterion Value"
  ) +
  annotate("text", x = best_aic, y = min(results$AIC),
           label = paste("Best Lag =", best_aic),
           vjust = 0, hjust = 1.2, color = "red", size = 2) +
  annotate("text", x = best_bic, y = min(results$BIC),
           label = paste("Best Lag =", best_aic),
           vjust = -4, hjust = 1, color = "steelblue", size = 2) +
  scale_x_continuous(limits = c(1, 12), breaks = 1:12) +
  theme_minimal()

# plot BIC
ggplot(results, aes(x = p, y = BIC)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 4, linetype = "dashed", color = "red") +
  labs(title = "BIC across AR(p) models", x = "Lag order (p)", y = "BIC") +
  theme_minimal()

# --------- AR(P) MODEL WITH ROLLING WINDOW ---------
# set seed for reproducibility
set.seed(4308)

# rolling window ar(p) forecast function
ar_rolling_window <- function(y, nprev = 66, p = 4, h = 1) {
  n <- length(y)
  preds   <- rep(NA_real_, nprev)
  actuals <- rep(NA_real_, nprev)
  aic_vals <- numeric(nprev)
  bic_vals <- numeric(nprev)
  
  for (i in 1:nprev) {
    # forecast origin index t (last index included in training)
    # i = 1 should forecast y[n - nprev + 1], so t = n - nprev + 1 - h
    t <- (n - nprev) + i - h
    if (t <= p) next  # not enough obs to fit AR(p)
    
    # training sample up to t
    y_train <- y[1:t]
    
    # build lag matrix for AR(p)
    lagged <- embed(y_train, p + 1)
    y_t <- lagged[, 1]
    X_t <- as.data.frame(lagged[, -1, drop = FALSE])
    colnames(X_t) <- paste0("lag", 1:p)
    
    train_df <- data.frame(y_t = y_t, X_t)
    model <- lm(y_t ~ ., data = train_df)
    
    # store ICs for this window (same across horizons by design)
    aic_vals[i] <- AIC(model)
    bic_vals[i] <- BIC(model)
    
    # recursive h-step ahead forecast:
    # keep lags in order [lag1=y_t, lag2=y_{t-1}, ..., lagp=y_{t-p+1}]
    lags_rev <- rev(tail(y_train, p))   # [y_t, y_{t-1}, ..., y_{t-p+1}]
    yhat_h <- NA_real_
    for (k in 1:h) {
      newdata <- as.data.frame(t(lags_rev))
      colnames(newdata) <- paste0("lag", 1:p)
      yhat_h <- as.numeric(predict(model, newdata = newdata))
      # update lags: new most-recent becomes yhat_h
      if (p > 1) {
        lags_rev <- c(yhat_h, lags_rev[1:(p-1)])
      } else {
        lags_rev <- c(yhat_h)
      }
    }
    
    preds[i]   <- yhat_h            # h-step prediction for target at t + h
    actuals[i] <- y[t + h]          # align actual at the same horizon
    
    cat("Rolling Window:", i, "/", nprev, "| p =", p, "| h =", h, "\n")
  }
  
  # metrics
  rmse_val <- sqrt(mean((preds - actuals)^2, na.rm = TRUE))
  mae_val  <- mean(abs(preds - actuals), na.rm = TRUE)
  
  list(
    preds = preds,
    actuals = actuals,
    rmse = rmse_val,
    mae  = mae_val,
    aic_mean = mean(aic_vals, na.rm = TRUE),
    bic_mean = mean(bic_vals, na.rm = TRUE),
    aic_series = aic_vals,
    bic_series = bic_vals
  )
}

# run the model
nprev <- 66
p <- 4

ar_roll <- ar_rolling_window(y = y_ts, nprev = nprev, p = p)

cat("Rolling-Window RMSE:", ar_roll$rmse, "\n")

# combine forecasts and actuals into one dataframe
plot_df <- data.frame(
  Index = seq_along(ar_roll$actuals),
  Actual = ar_roll$actuals,
  Forecast = ar_roll$preds
)

# convert to long format for ggplot
plot_df_long <- plot_df %>%
  tidyr::pivot_longer(cols = c("Actual", "Forecast"),
                      names_to = "Type",
                      values_to = "Value")

# plot
ggplot(plot_df_long, aes(x = Index, y = Value, color = Type)) +
  geom_line() +
  scale_color_manual(values = c("Actual" = "black", "Forecast" = "red")) +
  labs(
    title = paste0("Rolling AR(", p, ") Forecasts vs Actual Inflation"),
    subtitle = "1-Month Ahead Rolling Forecasts",
    x = "Time (last 66 months)",
    y = "Transformed Inflation Rate (CPI_t)",
    color = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )


# --------- LOOP ACROSS FORECASTS HORIZONS ---------
# set parameters
nprev <- 100
forecast_horizons <- c(1, 3, 6, 12)
lags <- 1:12

# directory to save outputs (optional)
dir.create("../data/ar(p)_results", showWarnings = FALSE)

# main loop
for (h in forecast_horizons) {
  for (p in lags) {
    
    cat("\nRunning AR(", p, ") with forecast horizon =", h, "...\n")
    
    ar_roll <- ar_rolling_window(y = y_ts, nprev = nprev, p = p, h = h)
    
    # name convention: horizon_1_ar1.rds, horizon_3_ar4.rds, etc.
    file_name <- sprintf("../data/ar(p)_results/horizon_%d_ar%d.rds", h, p)
    saveRDS(ar_roll, file = file_name)
    
    cat("Saved:", file_name, "\n")
  }
}


# --------- AGGREGATE RESULTS ---------
# define horizons and storage list
forecast_horizons <- c(1, 3, 6, 12)
results_list <- list()

# loop through each horizon
for (h in forecast_horizons) {
  # list all .rds files for this horizon
  files <- list.files("../data/ar(p)_results", 
                      pattern = paste0("horizon_", h, "_ar"), 
                      full.names = TRUE)
  
  horizon_df <- data.frame()  # empty df for that horizon
  
  for (file in files) {
    ar_roll <- readRDS(file)
    
    # extract lag number (from filename)
    p <- as.numeric(str_extract(file, "(?<=_ar)\\d+"))
    
    # append metrics
    horizon_df <- rbind(
      horizon_df,
      data.frame(
        Model = paste0("AR(", p, ")"),
        RMSE = ar_roll$rmse,
        MAE = ar_roll$mae,
        AIC_mean = ar_roll$aic_mean,
        BIC_mean = ar_roll$bic_mean
      )
    )
  }
  
  # arrange by lag order
  horizon_df <- horizon_df %>% arrange(as.numeric(str_extract(Model, "\\d+")))
  
  # store in list
  results_list[[paste0("H", h)]] <- horizon_df
}

results_list$H1
results_list$H3
results_list$H6
results_list$H12

# save the entire results list
saveRDS(results_list, file = "../data/ar(p)_results/summary_results_list.rds")

# find best ar(p) model across the different forecast horizonz
summary_best <- data.frame()

for (h in c(1, 3, 6, 12)) {
  df <- results_list[[paste0("H", h)]]
  
  summary_best <- rbind(
    summary_best,
    data.frame(
      Horizon = h,
      Best_RMSE_Model = df$Model[which.min(df$RMSE)],
      Best_MAE_Model  = df$Model[which.min(df$MAE)],
      Best_AIC_Model  = df$Model[which.min(df$AIC_mean)],
      Best_BIC_Model  = df$Model[which.min(df$BIC_mean)]
    )
  )
}

summary_best

# save results
saveRDS(summary_best, file = "../data/ar(p)_results/summary_best_model.rds")
