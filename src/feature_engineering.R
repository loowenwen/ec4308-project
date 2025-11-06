library(stats)
library(zoo)
library(tidyverse)
library(purrr)

# ---------------------- 1. Helper Functions -------------------------

# Create lag variables
create_lags <- function(df, n_lags, dates) {
  lagged_list <- list(df)  # Start with the original (unlagged) variables
  colnames(lagged_list[[1]]) <- colnames(df)  # ensure names are preserved
  
  # Add lagged versions
  for (l in 1:n_lags) {
    lagged <- dplyr::lag(df, l)
    colnames(lagged) <- paste0(colnames(df), "_lag", l)
    lagged_list[[l + 1]] <- lagged
  }
  
  lagged_df <- do.call(cbind, lagged_list)
  lagged_df <- cbind(date = dates, lagged_df)
  return(lagged_df)
}

# MARX transformation
create_marx <- function(df, max_lag, dates) {
  stopifnot(nrow(df) == length(dates))
  n <- nrow(df); k <- ncol(df)
  marx_matrix <- matrix(NA, nrow = n - max_lag, ncol = k * max_lag)
  
  for (p in 1:max_lag) {
    lagged_df <- dplyr::lag(df, p)
    marx_matrix[, ((p - 1) * k + 1):(p * k)] <- tail(as.matrix(lagged_df), n - max_lag)
  }
  
  marx_df <- as.data.frame(marx_matrix)
  colnames(marx_df) <- paste0(rep(colnames(df), each = max_lag), "_MARX", 1:max_lag)
  
  tibble(sasdate = tail(dates, n - max_lag), marx_df)
}

# MAF transformation
create_maf <- function(df, n_lags = 12, n_pcs = 1, dates) {
  stopifnot(nrow(df) == length(dates))
  n <- nrow(df)
  maf_list <- vector("list", ncol(df))
  
  for (k_idx in seq_along(df)) {
    x <- as.matrix(df[, k_idx])
    lagged_mat <- embed(x, n_lags)
    if (nrow(lagged_mat) < n_pcs) next
    
    pca_res <- prcomp(lagged_mat, scale. = TRUE)
    pcs <- pca_res$x[, 1:n_pcs, drop = FALSE]
    maf_list[[k_idx]] <- pcs
  }
  
  maf_matrix <- do.call(cbind, maf_list)
  colnames(maf_matrix) <- paste0(rep(colnames(df), each = n_pcs), "_MAF_PC", 1:n_pcs)
  
  tibble(
    sasdate = tail(dates, nrow(maf_matrix)),
    as.data.frame(maf_matrix)
  )
}

# ---------------------- 2. Combine train/test -------------------------

X_full <- rbind(X_train, X_test)
y_full <- rbind(y_train, y_test)

# ---------------------- 3. Prepare transformations -------------------------

p_y <- 4
p_m <- 4

X_t <- X_full %>% select(-sasdate)
X_t_raw_withdate <- fred_clean %>% select(-CPIAUCSL)
dates_Xt <- X_full$sasdate
dates_Xtraw <- X_t_raw_withdate$sasdate
dates_yt <- y_full$sasdate
X_t_raw <- X_t_raw_withdate %>% select(-sasdate)
y_t <- y_full

# ---- Lagged y and X ----
y_lags <- create_lags(y_t %>% select(-sasdate), p_y, dates_yt)
X_t_lags <- create_lags(X_t, p_m, dates_Xt)

# ---- MARX & MAF ----
marx_data <- create_marx(X_t, max_lag = 12, dates_Xt)
maf_data <- create_maf(X_t, n_lags = 12, n_pcs = 2, dates_Xt)

# ---------------------- 4. Save static components -------------------------

saveRDS(
  list(
    X_t = X_t,
    X_t_raw = X_t_raw,
    y_t = y_t,
    y_lags = y_lags,
    X_t_lags = X_t_lags,
    marx_data = marx_data,
    maf_data = maf_data,
    dates_Xt = dates_Xt,
    dates_Xtraw = dates_Xtraw
  ),
  file = "base_features_static.rds"
)

message("✅ Static feature engineering complete. Ready for rolling PCA.")




# Rolling-window PCA function
rolling_pca_lags <- function(X_full, dates_full, train_end_idx, max_lags = 4, var_threshold = 0.8, min_train = 50) {
  # Subset training window
  X_train <- X_full[1:train_end_idx, , drop = FALSE]
  if (nrow(X_train) < min_train) return(NULL)
  
  # PCA on training data
  pca_res <- prcomp(X_train, scale. = TRUE)
  cum_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_factors <- max(1, min(which(cum_var >= var_threshold)))
  
  # Factor scores for training period
  F_train <- predict(pca_res)[, 1:n_factors, drop = FALSE]
  F_train <- as.data.frame(F_train)
  colnames(F_train) <- paste0("F", 1:n_factors)
  F_train <- cbind(date = dates_full[1:train_end_idx], F_train)
  
  # Add lags including original factors
  lagged_factors <- create_lags(F_train[, -1, drop = FALSE], max_lags, dates_full[1:train_end_idx])
  return(lagged_factors)
}


#rolling build Z matrice
rolling_Z_builder <- function(X_t, X_t_raw, y_t, marx_data, 
                              dates_Xt, dates_Xtraw,
                              Z_names = c("Z_F_stationary", "Z_F_raw",
                                          "Z_X", "Z_X_MARX",
                                          "Z_Fstationary_X_MARX_","Z_Fraw_X_MARX_"),
                              train_window = 100, max_lags = 4, var_threshold = 0.8, min_train = 50) {
  
  n <- nrow(X_t)
  all_Z <- vector("list", n - train_window) # store Z matrices for each forecast date
  names(all_Z) <- paste0("t+", 1:(n - train_window))
  
  # Helper: align by date
  align_by_date <- function(...) {
    dfs <- list(...)
    min_rows <- min(sapply(dfs, nrow))
    dfs_trimmed <- lapply(dfs, function(x) x[(nrow(x)-min_rows+1):nrow(x), , drop=FALSE])
    Reduce(function(x,y) cbind(x, y[ , setdiff(names(y), "sasdate")]), dfs_trimmed)
  }
  
  for (t_idx in (train_window + 1):n) {
    train_start <- t_idx - train_window
    train_end <- t_idx - 1
    
    # Subset training window
    X_train <- X_t[train_start:train_end, , drop = FALSE]
    X_raw_train <- X_t_raw[train_start:train_end, , drop = FALSE]
    y_train <- y_t[train_start:train_end, , drop = FALSE]
    marx_train <- marx_data[train_start:train_end, , drop = FALSE]
    dates_train_Xt <- dates_Xt[train_start:train_end]
    dates_train_Xtraw <- dates_Xtraw[train_start:train_end]
    
    Z_step <- list()
    
    for (Z_name in Z_names) {
      if (Z_name == "Z_F_stationary") {
        F_lags <- rolling_pca_lags(X_train, dates_train_Xt, nrow(X_train), max_lags, var_threshold, min_train)
        Z_step[[Z_name]] <- align_by_date(F_lags, y_lags = create_lags(y_train %>% select(-sasdate), max_lags, dates_train_Xt))
        
      } else if (Z_name == "Z_F_raw") {
        F_lags <- rolling_pca_lags(X_raw_train, dates_train_Xtraw, nrow(X_raw_train), max_lags, var_threshold, min_train)
        Z_step[[Z_name]] <- align_by_date(F_lags, y_lags = create_lags(y_train %>% select(-sasdate), max_lags, dates_train_Xtraw))
        
      } else if (Z_name == "Z_X") {
        Z_step[[Z_name]] <- align_by_date(create_lags(X_train, max_lags, dates_train_Xt),
                                          create_lags(y_train %>% select(-sasdate), max_lags, dates_train_Xt))
        
      } else if (Z_name == "Z_X_MARX") {
        Z_step[[Z_name]] <- align_by_date(create_lags(X_train, max_lags, dates_train_Xt),
                                          marx_train,
                                          create_lags(y_train %>% select(-sasdate), max_lags, dates_train_Xt))
        
      } else if (Z_name == "Z_Fstationary_X_MARX_") {
        F_lags <- rolling_pca_lags(X_train, dates_train_Xt, nrow(X_train), max_lags, var_threshold, min_train)
        Z_step[[Z_name]] <- align_by_date(F_lags,
                                          create_lags(X_train, max_lags, dates_train_Xt),
                                          marx_train,
                                          create_lags(y_train %>% select(-sasdate), max_lags, dates_train_Xt))
        
      } else if (Z_name == "Z_Fraw_X_MARX_") {
        F_lags <- rolling_pca_lags(X_raw_train, dates_train_Xtraw, nrow(X_raw_train), max_lags, var_threshold, min_train)
        Z_step[[Z_name]] <- align_by_date(F_lags,
                                          create_lags(X_train, max_lags, dates_train_Xt),
                                          marx_train,
                                          create_lags(y_train %>% select(-sasdate), max_lags, dates_train_Xt))
      }
    }
    
    all_Z[[paste0("t+", t_idx)]] <- Z_step
  }
  
  return(all_Z)
}


