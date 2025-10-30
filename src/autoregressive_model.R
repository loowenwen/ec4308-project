## AR MODELS WITH Y LAGS ##

# load cleaned FRED data
fred_data_clean <- read.csv("../data/fred_train.csv")

fred_data_clean$sasdate <- as.Date(fred_data_clean$sasdate)

sample_start <- as.Date("2000-01-01")
sample_end <- as.Date("2019-12-01")

sample_period <- subset(
  fred_data_clean,
  sasdate >= sample_start & sasdate <= sample_end
)

Y <- sample_period$CPIAUCSL

if (any(is.na(Y))) {
  stop("Sample contains missing values in CPIAUCSL.")
}

RMSE <- function(pred, truth) {
  sqrt(mean((truth - pred)^2))
}

nprev <- 70
max_lag <- 12

if (length(Y) <= (max_lag + nprev)) {
  stop("Not enough observations for the requested lag structure and forecast window.")
}

lagged <- embed(Y, max_lag + 1)
targets <- lagged[, 1]
lag_matrix <- lagged[, -1, drop = FALSE]
colnames(lag_matrix) <- paste0("lag", seq_len(max_lag))

total_rows <- nrow(lagged)
test_rows <- seq(total_rows - nprev + 1, total_rows)
train_rows <- seq_len(total_rows - nprev)

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

  rmse_val <- RMSE(preds, test_df$y)
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
  "Sample period:", sample_start, "to", sample_end, "\n",
  "Observations used:", length(Y), "\n",
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
