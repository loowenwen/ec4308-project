# ---------- LOAD & SETUP ----------
library(dplyr)
library(lubridate)
library(readr)
library(tidyverse)
library(randomForest)
set.seed(4308)

# ---------- PART 1 ----------

# ---------- HELPER: LEAD TARGET ----------
make_horizon_target <- function(y_df, h) {
  y_df %>%
    arrange(sasdate) %>%
    mutate(target = dplyr::lead(CPI_t, h)) %>%
    select(sasdate, target)
}

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



# ---------- PART 2 ----------

# ---------- HELPER: LEAD TARGET ----------
make_horizon_target <- function(y_df, h) {
  y_df %>%
    arrange(sasdate) %>%
    mutate(target = dplyr::lead(CPI_t, h)) %>%
    select(sasdate, target)
}

# ---------- HELPER: CREATE LAGS ----------
create_lags <- function(df, n_lags, dates) {
  lagged_list <- list(df)
  for (l in 1:n_lags) {
    lagged <- dplyr::lag(df, l)
    colnames(lagged) <- paste0(colnames(df), "_lag", l)
    lagged_list[[l + 1]] <- lagged
  }
  lagged_df <- do.call(cbind, lagged_list)
  cbind(sasdate = dates, lagged_df)
}

# ---------- HELPER: ALIGN BY MOST RECENT OBSERVATIONS ----------
align_by_date_join <- function(...) {
  dfs <- list(...)
  out <- dfs[[1]]
  stopifnot("sasdate" %in% names(out))
  for (k in 2:length(dfs)) {
    stopifnot("sasdate" %in% names(dfs[[k]]))
    out <- dplyr::left_join(out, dfs[[k]], by = "sasdate")
  }
  out
}

# ---------- HELPER: ROLLING PCA LAGS ----------
rolling_pca_lags <- function(X_full, dates_full, train_end_idx,
                             max_lags = 4, var_threshold = 0.8, min_train = 50) {
  # subset data up to train_end_idx
  X_train <- X_full[1:train_end_idx, , drop = FALSE]
  
  # keep only numeric columns (avoid colMeans error)
  X_train <- X_train %>% select(where(is.numeric)) %>% as.data.frame()
  
  # check minimum length
  if (nrow(X_train) < min_train) {
    warning("Too few rows for PCA in this window.")
    return(NULL)
  }
  
  # PCA
  pca_res <- prcomp(X_train, scale. = TRUE)
  
  # determine number of factors to keep
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_factors <- max(1, min(which(cum_var >= var_threshold)))
  
  # extract factor scores for training period
  F_train <- predict(pca_res)[, 1:n_factors, drop = FALSE] %>% as.data.frame()
  colnames(F_train) <- paste0("F", 1:n_factors)
  F_train <- cbind(sasdate = dates_full[1:train_end_idx], F_train)
  
  # create lagged factors
  lagged_factors <- create_lags(F_train[, -1, drop = FALSE],
                                max_lags,
                                dates_full[1:train_end_idx])
  
  # return final lagged factor dataframe
  return(tail(lagged_factors, nrow(X_train) - max_lags))
}

# ---------- HELPER: ROLLING RANDOM FOREST WITH PCA ----------
rf_rolling_window_pca <- function(Z_naked, X_stationary, y, dates_X,
                                  nprev = 100, h = 1, ntree = 1000, mtry = NULL,
                                  max_lags = 4, var_threshold = 0.8) {
  
  y_h <- make_horizon_target(y, h)
  
  df <- Z_naked %>%
    left_join(y_h, by = "sasdate") %>%
    arrange(sasdate) %>%
    drop_na()
  
  n <- nrow(df)
  train_end <- n - nprev
  
  # storage
  preds_raw <- preds_stat <- actuals <- dates <- numeric(nprev)
  
  cat("Rolling RF (Stat PCA) | Horizon =", h, "| OOS =", nprev, "\n")
  
  for (i in 1:nprev) {
    train_idx <- 1:(train_end + i - 1)
    test_idx  <- train_end + i
    
    # build factor lags on TRAIN ONLY
    F_stat_lags <- rolling_pca_lags(X_stationary[train_idx, ], dates_X[train_idx],
                                    train_end_idx = length(train_idx),
                                    max_lags, var_threshold)
    
    # if PCA failed early in sample (too short), skip this step
    if (is.null(F_stat_lags)) {
      preds_stat[i] <- NA_real_
      actuals[i]    <- df$target[test_idx]
      dates[i]      <- df$sasdate[test_idx]
      next}
    
    # build Z up to test point (1:test_idx) so "last row" is the test observation
    Z_full_stat <- align_by_date_join(Z_naked[1:test_idx, c("sasdate", setdiff(names(Z_naked), "sasdate"))], F_stat_lags)
    
    # keep only complete rows (avoid NA in predictors)
    Z_full_stat <- Z_full_stat[complete.cases(Z_full_stat), ]
    
    # after complete.cases(), the last row is the test row; split train/test
    Z_train_stat <- Z_full_stat[-nrow(Z_full_stat), , drop = FALSE]
    Z_test_stat  <- Z_full_stat[nrow(Z_full_stat), , drop = FALSE]
    
    # features (drop sasdate)
    feature_cols_stat <- setdiff(names(Z_train_stat), "sasdate")
    train_x_stat      <- Z_train_stat[, feature_cols_stat, drop = FALSE]
    test_x_stat       <- Z_test_stat[, feature_cols_stat, drop = FALSE]
    
    # align y to X (use tail so lengths match exactly after shrinkage)
    train_y_stat <- tail(df$target[1:(test_idx - 1)], nrow(train_x_stat))
    test_y       <- df$target[test_idx]
    
    # random forest (stationary factors)
    rf_stat <- randomForest(
      x = train_x_stat, y = train_y_stat,
      ntree = ntree,
      mtry = ifelse(is.null(mtry), floor(sqrt(ncol(train_x_stat))), mtry)
    )
    preds_stat[i] <- predict(rf_stat, test_x_stat)
    
    actuals[i] <- test_y
    dates[i]   <- df$sasdate[test_idx]
    
    if (i %% 10 == 0 || i == 1) {
      cat(sprintf("[%3d/%3d] %s | Stat=%.4f | Actual=%.4f\n",
                  i, nprev, as.character(dates[i]),
                  preds_stat[i], actuals[i]))
    }
  }
  
  # compute metrics
  res_stat <- tibble(sasdate = as.Date(dates, origin = "1970-01-01"),
                     yhat = preds_stat, y = actuals)
  
  metrics_stat <- tibble(H = h,
                         RMSE = sqrt(mean((res_stat$yhat - res_stat$y)^2, na.rm = TRUE)),
                         MAE  = mean(abs(res_stat$yhat - res_stat$y), na.rm = TRUE))
  
  cat(sprintf("RMSE_stat = %.4f\n", metrics_stat$RMSE))
  
  list(results = res_stat, metrics = metrics_stat)
}


# ---------- USING FEATURE ENGINEERED DATA ----------
# load precomputed Z matrices and predictors
Z_list <- readRDS("all_Z_matrices.rds")

# inspect what’s inside
names(Z_list)

# unpack individual z-matrices
# note: decided to not proceed with MAF features
Z_F_naked <- Z_list$Z_F_naked
Z_F_Level_naked <- Z_list$Z_F_Level_naked
Z_X <- Z_list$Z_X
Z_Ht <- Z_list$Z_Ht
Z_X_MARX <- Z_list$Z_X_MARX
Z_F_X_MARX_naked <- Z_list$Z_F_X_MARX_naked
Z_F_X_MARX_Level_naked <- Z_list$Z_F_X_MARX_Level_naked

# load raw & stationary predictors for PCA
X_stationary <- read_csv("../data/predictors_transformed.csv") %>% select(-sasdate)
dates_X <- read_csv("../data/predictors_transformed.csv")$sasdate

# target
y <- read_csv("../data/cpi_target_full.csv") %>%
  select(sasdate, CPI_t)



# ---------- Z_F_naked FEATURE ----------
Z_naked <- Z_F_naked

# ensure every Z has a sasdate column aligned to y
if (!"sasdate" %in% names(Z_naked)) {
  # length difference between y and Z_naked
  n_diff <- length(y$sasdate) - nrow(Z_naked)
  
  # assign the last nrow(Z_naked) dates, ending at 2025-06
  Z_naked <- cbind(
    sasdate = tail(y$sasdate, nrow(Z_naked)),
    Z_naked
  )
}

# run for one horizon
set.seed(4308)
rf_results_h1 <- rf_rolling_window_pca(
  Z_naked = Z_naked,
  X_stationary = X_stationary,
  y = y,
  dates_X = dates_X,
  nprev = 100,
  h = 1
)

# save results
saveRDS(rf_results_h1, "../data/rf_results/Z_Fstationary_h1.rds")

# run for remaining horizon using a loop
horizons <- c(3, 6, 12)

set.seed(4308)
for (h in horizons) {
  cat("\nRunning RF with PCA (Stationary) | Horizon =", h, "\n")
  
  rf_results <- rf_rolling_window_pca(
    Z_naked = Z_naked,
    X_stationary = X_stationary,
    y = y,
    dates_X = dates_X,
    nprev = 100,
    h = h
  )
  
  # save each result with horizon in filename
  save_path <- paste0("../data/rf_results/Z_Fstationary_h", h, ".rds")
  saveRDS(rf_results, save_path)
  
  cat("Saved results to:", save_path, "\n")
}


# ---------- LOOP THROUGH ALL Z'S THAT DOES FACTORS ON STATIONARY X's ----------
Z_list_to_run <- list(
  Z_Level_F        = Z_list$Z_F_Level_naked,
  Z_X_MARX_F       = Z_list$Z_F_X_MARX_naked,
  Z_X_MARX_Level_F = Z_list$Z_F_X_MARX_Level_naked
)

# loop settings
horizons <- c(1, 3, 6, 12)
set.seed(4308)

# main loop
for (Z_name in names(Z_list_to_run)) {
  cat("\n========================================\n")
  cat("Running feature set:", Z_name, "\n")
  cat("========================================\n")
  
  Z_naked <- Z_list_to_run[[Z_name]]
  
  # ensure every Z has a sasdate column aligned to y
  if (!"sasdate" %in% names(Z_naked)) {
    n_diff <- length(y$sasdate) - nrow(Z_naked)
    Z_naked <- cbind(
      sasdate = tail(y$sasdate, nrow(Z_naked)),
      Z_naked
    )
  }
  
  # run across all horizons
  for (h in horizons) {
    cat("\n--- Horizon", h, "months ahead ---\n")
    
    rf_results <- rf_rolling_window_pca(
      Z_naked = Z_naked,
      X_stationary = X_stationary,
      y = y,
      dates_X = dates_X,
      nprev = 100,
      h = h
    )
    
    # clean filename
    save_path <- paste0("../data/rf_results/", Z_name, "stationary_h", h, ".rds")
    saveRDS(rf_results, save_path)
    
    cat("Saved results to:", save_path, "\n")
  }
}



# ---------- LOOP THROUGH ALL Z'S THAT DOES NOT HAVE FACTORS ON STATIONARY X's ----------
# load the y values
cpi_target_full <- read_csv("../data/cpi_target_full.csv")
y <- cpi_target_full[, -2]  # assuming column 2 is dropped

Z_list_to_run <- list(
  Z_X        = Z_list$Z_X, 
  Z_Ht       = Z_list$Z_Ht, 
  Z_X_MARX   = Z_list$Z_X_MARX
)

# loop settings
horizons <- c(1, 3, 6, 12)
set.seed(4308)

# main loop
for (Z_name in names(Z_list_to_run)) {
  cat("\n========================================\n")
  cat("Running Feature Set:", Z_name, "\n")
  cat("========================================\n")
  
  Z <- Z_list_to_run[[Z_name]]
  
  # run across all horizons
  for (h in horizons) {
    cat("\nStarting", h, "month(s) ahead forecast at", format(Sys.time(), "%H:%M:%S"), "\n")
    
    rf_results <- rf_rolling_window(Z, y, nprev = 100, h = h)
    
    # clean filename
    save_path <- paste0("../data/rf_results/", Z_name, "_h", h, ".rds")
    saveRDS(rf_results, save_path)
    
    cat("Saved results to:", save_path, "\n")
    cat("\nCompleted", h, "month(s) ahead forecast at", format(Sys.time(), "%H:%M:%S"), "\n")
  }
}
