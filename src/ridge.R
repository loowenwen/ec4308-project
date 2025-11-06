library(dplyr)
library(glmnet)
library(doParallel)
library(foreach)


static_features <- readRDS("base_features_static.rds")
X_t         <- static_features$X_t          # transformed (stationary)
X_t_raw     <- static_features$X_t_raw      # raw levels
y_t         <- static_features$y_t          # CPI_t
y_lags      <- static_features$y_lags       # lagged CPI
X_t_lags    <- static_features$X_t_lags     # lagged predictors
marx_data   <- static_features$marx_data    # extra variables

## --------------------------------------------------------------
## 2. Helper functions
## --------------------------------------------------------------
MSE  <- function(pred, truth) mean((pred - truth)^2, na.rm = TRUE)
RMSE <- function(pred, truth) sqrt(MSE(pred, truth))
MAE  <- function(pred, truth) mean(abs(pred - truth), na.rm = TRUE)

rolling_pca_lags <- function(X_full, dates_full, train_end_idx, p_f = 4) {
  X_train <- X_full[1:train_end_idx, , drop = FALSE]
  if (nrow(X_train) < 50) return(NULL)
  
  # PCA on training window
  pca_res <- prcomp(X_train, scale. = TRUE)
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_factors <- max(1, min(which(cum_var >= 0.8)))
  
  # Factor scores for training period
  F_train <- predict(pca_res)[, 1:n_factors, drop = FALSE]
  F_train <- as.data.frame(F_train)
  colnames(F_train) <- paste0("F", 1:n_factors)
  F_train <- cbind(sasdate = dates_full[1:train_end_idx], F_train)
  
  # Add current + lags
  create_lags(F_train[, -1], p_f, dates_full[1:train_end_idx])
}

clean_Z_matrix <- function(Z_mat, target_col = "CPI_t") {
  if (!target_col %in% colnames(Z_mat)) stop(paste("Missing target:", target_col))
  Z_mat <- Z_mat[!is.na(Z_mat[[target_col]]), ]
  pred_cols <- setdiff(colnames(Z_mat), c("sasdate", target_col))
  Z_mat <- Z_mat[!apply(Z_mat[, pred_cols, drop = FALSE], 1, anyNA), ]
  Z_mat <- Z_mat[, !duplicated(colnames(Z_mat))]
  rownames(Z_mat) <- NULL
  Z_mat
}

create_lags <- function(df, p, dates) {
  
  stopifnot(is.data.frame(df), length(dates) >= nrow(df))
  
  ## Keep only rows that have a date
  valid <- 1:nrow(df)
  if (length(dates) > nrow(df)) {
    dates <- dates[1:nrow(df)]  # in case dates is longer
  }
  
  if (p == 0) {
    return(tibble(sasdate = dates, df))
  }
  
  ## Create lagged versions
  lagged_list <- lapply(1:p, function(l) dplyr::lag(df, l))
  lagged_df   <- do.call(cbind, lagged_list)
  colnames(lagged_df) <- paste0(rep(colnames(df), each = p), "_L", rep(1:p, times = ncol(df)))
  
  ## Combine current + lags
  full_df <- cbind(df, lagged_df)
  
  ## Valid rows: after p lags, no NA in lagged columns
  valid_rows <- (p + 1):nrow(df)
  
  tibble(
    sasdate = dates[valid_rows],
    as.data.frame(full_df[valid_rows, , drop = FALSE])
  )
}

## --------------------------------------------------------------
## 3. Rolling-window PCA helper
## --------------------------------------------------------------
roll_pca_factors <- function(X_full, dates_full, train_end_idx, p_f = 4) {
  X_train <- X_full[1:train_end_idx, , drop = FALSE]
  if (nrow(X_train) < 50) return(NULL)  # safety
  
  pca_res <- prcomp(X_train, scale. = TRUE)
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_factors <- min(which(cum_var >= 0.8))
  
  F_all <- predict(pca_res, newdata = X_full)[, 1:n_factors, drop = FALSE]
  F_all <- as.data.frame(F_all)
  colnames(F_all) <- paste0("F", 1:n_factors)
  
  F_all <- cbind(sasdate = dates_full, F_all)
  F_lags <- create_lags(F_all[, -1], p_f, dates_full)
  F_lags
}

## --------------------------------------------------------------
## 4. Build raw data pieces
## --------------------------------------------------------------
build_full_Z <- function(Z_type, p_f = 4) {
  has_date <- function(df) "sasdate" %in% colnames(df)
  inputs <- list(y_lags = y_lags, X_t_lags = X_t_lags,
                 marx = marx_data, y_t = y_t)
  if (grepl("^Z_F", Z_type)) {
    inputs$X_t <- X_t; inputs$X_t_raw <- X_t_raw
  }
  
  date_cols <- lapply(inputs, function(df) if (has_date(df)) df$sasdate else NULL)
  date_cols <- date_cols[!sapply(date_cols, is.null)]
  if (length(date_cols) == 0) stop("No sasdate found")
  common_dates <- Reduce(intersect, date_cols)
  
  subset_to_common <- function(df) {
    if (has_date(df)) {
      df %>% filter(sasdate %in% common_dates) %>% arrange(sasdate)
    } else {
      df[seq_along(common_dates), , drop = FALSE]
    }
  }
  
  y_lags_win <- subset_to_common(y_lags)
  X_lags_win <- subset_to_common(X_t_lags)
  marx_win   <- subset_to_common(marx_data)
  y_win      <- subset_to_common(y_t)
  X_win <- X_win_raw <- NULL
  if (grepl("^Z_F", Z_type)) {
    X_win     <- subset_to_common(X_t) %>% select(-any_of("sasdate"))
    X_win_raw <- subset_to_common(X_t_raw) %>% select(-any_of("sasdate"))
  }
  dates_win <- y_lags_win$sasdate
  
  list(
    y_lags_win = y_lags_win,
    X_lags_win = X_lags_win,
    marx_win   = marx_win,
    y_win      = y_win,
    X_win      = X_win,
    X_win_raw  = X_win_raw,
    dates_win  = dates_win,
    Z_type     = Z_type,
    p_f        = p_f
  )
}


## --------------------------------------------------------------
## 5. Accurate Ridge Forecast – Rolling PCA + Purged CV (FINAL)
## --------------------------------------------------------------
accurate_ridge_forecast <- function(Z_parts, nprev = 100, horizons = c(1,3,6,12)) {

  y_raw   <- Z_parts$y_win$CPI_t
  dates   <- Z_parts$dates_win
  n_total <- length(y_raw)
  
  ## raw predictor matrix (only for factor models)
  X_raw <- if (grepl("^Z_F", Z_parts$Z_type))
    if (grepl("F_stationary", Z_parts$Z_type)) Z_parts$X_win
  else Z_parts$X_win_raw
  else NULL
  
  ## base matrix = lagged X / MARX + CPI target
  base_mat <- switch(Z_parts$Z_type,
                     "Z_X" = left_join(Z_parts$X_lags_win, Z_parts$y_lags_win, by = "sasdate"),
                     "Z_X_MARX" = left_join(Z_parts$X_lags_win, Z_parts$marx_win, by = "sasdate") %>%
                       left_join(Z_parts$y_lags_win, by = "sasdate"),
                     "Z_F_stationary" = Z_parts$y_lags_win,
                     "Z_F_raw"        = Z_parts$y_lags_win,
                     "Z_Fstationary_X_MARX_" = left_join(Z_parts$X_lags_win, Z_parts$marx_win, by = "sasdate") %>%
                       left_join(Z_parts$y_lags_win, by = "sasdate"),
                     "Z_Fraw_X_MARX_" = left_join(Z_parts$X_lags_win, Z_parts$marx_win, by = "sasdate") %>%
                       left_join(Z_parts$y_lags_win, by = "sasdate"),
                     stop("Invalid Z_type")
  )
  base_mat <- clean_Z_matrix(base_mat, "CPI_t")
  
  min_train <- 120 + 12 + 10
  first_train_end <- min_train
  if (first_train_end + nprev > nrow(base_mat)) stop("Not enough data")
  
  ## --------------------------------------------------------------
  ## 2. Loop over horizons
  ## --------------------------------------------------------------
  results <- list()
  for (h in horizons) {
    preds   <- actuals <- rep(NA, nprev)
    
    for (i in 1:nprev) {
      train_end_idx <- first_train_end + i - 1
      train_start   <- max(1, train_end_idx - 119)
      
      ## ----------------------------------------------------------
      ## 2.1  Build the training window (base + rolling PCA)
      ## ----------------------------------------------------------
      Z_full <- base_mat %>% slice(1:train_end_idx)   # <-- defined first
      
      F_lagged <- NULL
      if (grepl("^Z_F", Z_parts$Z_type)) {
        F_lagged <- rolling_pca_lags(X_raw, dates, train_end_idx,
                                     p_f = Z_parts$p_f)
        if (is.null(F_lagged)) next
        Z_full <- left_join(Z_full, F_lagged, by = "sasdate")
      }
      
      Z_train <- clean_Z_matrix(Z_full, "CPI_t")
      if (nrow(Z_train) < 100) next
      
      X_train <- as.matrix(Z_train %>% select(-sasdate, -CPI_t))
      y_train <- Z_train$CPI_t
      

      # --- 3. Purged CV (h-step target) – NEVER skip, standardize = TRUE ---
      y_cv <- y_train
      if (sd(y_train) == 0 || length(unique(y_train)) <= 1) {
        # Add *tiny* jitter only for CV (does NOT affect final model)
        y_cv <- y_train + rnorm(length(y_train), 0, 1e-8)
      }
      
      cv_folds  <- 5
      fold_size <- max(1, floor(nrow(Z_train) / (cv_folds + 1)))
      lambdas   <- 10^seq(3, -3, length = 50)
      mse_folds <- matrix(NA, nrow = cv_folds, ncol = length(lambdas))
      
      for (f in 1:cv_folds) {
        val_end   <- nrow(Z_train) - f * fold_size
        val_start <- val_end - fold_size + 1
        if (val_start <= 1) break
        
        train_cv <- X_train[1:(val_start-1), , drop = FALSE]
        y_train_cv <- y_cv[1:(val_start-1)]
        val_X    <- X_train[val_start:val_end, , drop = FALSE]
        
        val_y_idx <- (train_start + val_start - 1) + (h - 1) + (0:(val_end - val_start))
        val_y    <- dplyr::lead(y_raw, h)[val_y_idx]
        if (length(val_y) != nrow(val_X)) next
        
        # Jitter val_y if constant
        if (sd(val_y) == 0) val_y <- val_y + rnorm(length(val_y), 0, 1e-8)
        
        for (j in seq_along(lambdas)) {
          fit_cv <- glmnet(train_cv, y_train_cv,
                           alpha = 0, lambda = lambdas[j],
                           standardize = TRUE)  # ← ALWAYS TRUE
          pred   <- predict(fit_cv, val_X, s = lambdas[j])
          mse_folds[f, j] <- mean((pred - val_y)^2, na.rm = TRUE)
        }
      }
      best_lambda <- lambdas[which.min(colMeans(mse_folds, na.rm = TRUE))]
      

      # --- 4. Final model – jitter only if y constant ---
      y_final <- y_train
      if (sd(y_train) == 0 || length(unique(y_train)) <= 1) {
        y_final <- y_train + rnorm(length(y_train), 0, 1e-8)
      }
      
      fit <- glmnet(X_train, y_final, alpha = 0, lambda = best_lambda, standardize = TRUE)
      
      ## ----------------------------------------------------------
      ## 2.4  Build the test row (exact same column layout)
      ## ----------------------------------------------------------
      test_idx <- train_end_idx + 1
      if (test_idx > nrow(base_mat)) break
      
      Z_test <- base_mat[test_idx, , drop = FALSE]
      
      if (!is.null(F_lagged)) {
        ## reuse the PCA rotation from the training window
        pca_res   <- prcomp(X_raw[1:train_end_idx, , drop = FALSE],
                            scale. = TRUE)
        n_factors <- ncol(F_lagged) - 1   # current + 4 lags → 5 cols per factor
        n_factors <- n_factors %/% 5
        
        ## current factor scores for the test point
        F_cur <- predict(pca_res,
                         newdata = X_raw[test_idx, , drop = FALSE])[, 1:n_factors,
                                                                    drop = FALSE]
        F_cur <- as.data.frame(F_cur)
        colnames(F_cur) <- paste0("F", 1:n_factors)
        
        ## reconstruct the full lagged row (current + 4 lags)
        F_test_full <- F_cur
        for (lg in 1:Z_parts$p_f) {
          lag_cols <- paste0("F", 1:n_factors, "_l", lg)
          lag_vals <- F_lagged %>%
            filter(sasdate == dates[test_idx - lg]) %>%
            select(all_of(lag_cols))
          if (nrow(lag_vals) == 0) {
            lag_vals <- matrix(NA, nrow = 1, ncol = n_factors)
          }
          colnames(lag_vals) <- lag_cols
          F_test_full <- cbind(F_test_full, lag_vals)
        }
        Z_test <- cbind(Z_test, F_test_full)
      }
      
      Z_test <- clean_Z_matrix(Z_test, "CPI_t")
      if (nrow(Z_test) == 0) next
      
      X_test <- as.matrix(Z_test %>% select(-sasdate, -CPI_t))
      
      ## column-count safety
      if (ncol(X_test) != ncol(X_train)) next
      
      preds[i]   <- drop(predict(fit, newx = X_test, s = best_lambda))
      actual_idx <- test_idx + h - 1
      if (actual_idx <= n_total) actuals[i] <- y_raw[actual_idx]
    }
    
    results[[paste0("h", h)]] <- list(
      pred = preds, real = actuals,
      rmse = RMSE(preds, actuals), mae = MAE(preds, actuals)
    )
  }
  results
}

## --------------------------------------------------------------
## 6. Parallel setup
## --------------------------------------------------------------
n_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
options(warn = -1)
on.exit(stopCluster(cl), add = TRUE)

## --------------------------------------------------------------
## 7. Parameters & main loop
## --------------------------------------------------------------
Z_names <- c("Z_F_stationary","Z_F_raw","Z_X","Z_X_MARX",
             "Z_Fstationary_X_MARX_","Z_Fraw_X_MARX_")
horizons <- c(1, 3, 6, 12)
nprev <- 100
ridge_path <- "accurate_rolling_pca_ridge_results.rds"

all_ridge <- if (file.exists(ridge_path)) readRDS(ridge_path) else list()
Z_to_do <- setdiff(Z_names, names(all_ridge))

for (Z_type in Z_to_do) {
  cat("\n=== Rolling PCA + Ridge for", Z_type, "===\n")
  Z_parts <- build_full_Z(Z_type, p_f = 4)
  n_total <- nrow(Z_parts$y_win)
  cat(" Aligned rows:", n_total, "\n")
  
  if (n_total < 120 + 12 + 10 + nprev) {
    cat(" Not enough data, skipping\n")
    next
  }
  
  res <- accurate_ridge_forecast(Z_parts, nprev = nprev, horizons = horizons)
  all_ridge[[Z_type]] <- res
  saveRDS(all_ridge, ridge_path)
  cat(" Saved. RMSEs:",
      paste(sapply(horizons, function(h)
        paste0("h", h, "=", round(res[[paste0("h", h)]]$rmse, 4))), collapse = " | "),
      "\n")
}

cat("\nAll done. Results in:", ridge_path, "\n")


# Collect RMSEs for each horizon
mse_summary <- do.call(rbind, lapply(names(ridge_results), function(z_name) {
  res_list <- ridge_results[[z_name]]
  do.call(rbind, lapply(names(res_list), function(h) {
    rmse_val <- res_list[[h]]$rmse
    data.frame(
      Z_name = z_name,
      Horizon = h,
      RMSE = rmse_val,
      MSE = rmse_val^2
    )
  }))
})) %>%
  as.data.frame()

# Rank within each horizon (1 = best / smallest MSE)
mse_ranked <- mse_summary %>%
  group_by(Horizon) %>%
  mutate(Rank = rank(MSE, ties.method = "first")) %>%
  arrange(Horizon, Rank) %>%
  ungroup()

# Print and save results
print(mse_ranked)

fred_raw$G
ggplot(fred_clean, aes(x = date, y = GDP)) +
  geom_line(color = "#2C3E50", linewidth = 1) +
  labs(
    title = "US GDP Over Time",
    x = "Year",
    y = "GDP (in billions or appropriate units)"
  )
