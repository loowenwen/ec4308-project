# ---------- LOAD & SETUP ----------
library(dplyr)
library(lubridate)
library(randomForest)
set.seed(4308)

# ---------- HELPER: LEAD TARGET ----------
make_horizon_target <- function(y_df, h) {
  y_df %>%
    arrange(sasdate) %>%
    mutate(target = dplyr::lead(CPI_t, h)) %>%
    select(sasdate, target)
}

make_horizon_target(y, 1)

# ---------- HELPER: ROLLING RANDOM FOREST ----------
rf_rolling_window <- function(X, y, nprev = 100, h = 1, ntree = 1000, mtry = NULL) {
  # align target for horizon h
  y_h = make_horizon_target(y, h)
  
  # merge predictors and target
  df = X %>%
    left_join(y_h, by = "sasdate") %>%
    arrange(sasdate) %>%
    drop_na()
  
  n = nrow(df)
  train_end = n - nprev
  
  feature_cols = setdiff(names(df), c("sasdate", "target"))
  
  preds <- actuals <- dates <- numeric(nprev)
  
  cat("Rolling Window RF | Horizon =", h, "| OOS =", nprev, "\n")
  
  for (i in 1:nprev) {
    train_idx <- 1:(train_end + i - 1)
    test_idx  <- train_end + i
    
    train_x <- df[train_idx, feature_cols]
    train_y <- df[train_idx, "target", drop = TRUE]
    test_x  <- df[test_idx, feature_cols, drop = FALSE]
    test_y  <- df[test_idx, "target", drop = TRUE]
    
    rf <- randomForest(
      x = train_x,
      y = train_y,
      ntree = ntree,
      mtry = ifelse(is.null(mtry), floor(sqrt(ncol(train_x))), mtry),
      importance = TRUE
    )
    
    yhat <- predict(rf, test_x)
    
    preds[i] <- yhat
    actuals[i] <- test_y
    dates[i] <- df$sasdate[test_idx]
    
    if (i %% 10 == 0 || i == 1) {
      cat(sprintf("[%3d/%3d] %s | yhat = %.4f | actual = %.4f\n",
                  i, nprev, as.character(df$sasdate[test_idx]), yhat, test_y))
    }
  }
  
  results <- tibble(
    sasdate = as.Date(dates, origin = "1970-01-01"),
    yhat = preds,
    y = actuals
  )
  
  rmse <- sqrt(mean((results$yhat - results$y)^2, na.rm = TRUE))
  mae  <- mean(abs(results$yhat - results$y), na.rm = TRUE)
  
  cat(sprintf("Final RMSE = %.4f | MAE = %.4f\n", rmse, mae))
  
  list(
    results = results,
    metrics = tibble(H = h, RMSE = rmse, MAE = mae)
  )
}

# ---------- USING PREDICTORS_TRANSFORMED DATA ----------
# load the X values
predictors_transformed <- read_csv("../data/predictors_transformed.csv")
X <- predictors_transformed

# load the y values
cpi_target_full <- read_csv("../data/cpi_target_full.csv")
y <- cpi_target_full[, -2]  # assuming column 2 is dropped

# run for 1-month ahead forecast
rf_h1 <- rf_rolling_window(X, y, nprev = 100, h = 1)
saveRDS(rf_h1, "../data/rf_results/predictors_transformed_h1.rds")

# run for other horizons
rf_h3 <- rf_rolling_window(X, y, nprev = 100, h = 3)
saveRDS(rf_h3, "../data/rf_results/predictors_transformed_h3.rds")
rf_h6 <- rf_rolling_window(X, y, nprev = 100, h = 6)
saveRDS(rf_h6, "../data/rf_results/predictors_transformed_h6.rds")
rf_h12 <- rf_rolling_window(X, y, nprev = 100, h = 12)
saveRDS(rf_h12, "../data/rf_results/predictors_transformed_h12.rds")


# ---------- USING RAW DATA (NON-STATIONARY) ----------
# load the X values
# previously split on Jan 2020 
fred_train <- read_csv("../data/fred_train.csv")
fred_test <- read_csv("../data/fred_test.csv")
predictors_raw <- bind_rows(list(fred_train, fred_test))
X <- predictors_raw

# load the y values
cpi_target_full <- read_csv("../data/cpi_target_full.csv")
y <- cpi_target_full[, -2]  # assuming column 2 is dropped

# run for 1-month ahead forecast
rf_h1 <- rf_rolling_window(X, y, nprev = 100, h = 1)
saveRDS(rf_h1, "../data/rf_results/predictors_raw_h1.rds")
cat("\nCompleted 1-month ahead forecast at", format(Sys.time(), "%H:%M:%S"), "\n")

# run for other horizons
rf_h3 <- rf_rolling_window(X, y, nprev = 100, h = 3)
saveRDS(rf_h3, "../data/rf_results/predictors_raw_h3.rds")
cat("\nCompleted 3-month ahead forecast at", format(Sys.time(), "%H:%M:%S"), "\n")
rf_h6 <- rf_rolling_window(X, y, nprev = 100, h = 6)
saveRDS(rf_h6, "../data/rf_results/predictors_raw_h6.rds")
cat("\nCompleted 6-month ahead forecast at", format(Sys.time(), "%H:%M:%S"), "\n")
rf_h12 <- rf_rolling_window(X, y, nprev = 100, h = 12)
saveRDS(rf_h12, "../data/rf_results/predictors_raw_h12.rds")  
cat("\nCompleted 12-month ahead forecast at", format(Sys.time(), "%H:%M:%S"), "\n")