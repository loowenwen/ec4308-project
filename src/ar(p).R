# load libraries
library(forecast)
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
nprev <- 66 
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
  "Observations used:", length(y), "\n",
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