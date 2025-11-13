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
  lagged_df <- cbind(sasdate = dates, lagged_df)
  return(lagged_df)
}

#align by date function
align_by_date<- function(...) {
  dfs <- list(...)
  min_rows <- min(sapply(dfs, nrow))
  dfs_trimmed <- lapply(dfs, function(x) x[(nrow(x)-min_rows+1):nrow(x), , drop=FALSE])
  Reduce(function(x,y) cbind(x, y[ , setdiff(names(y), "sasdate")]), dfs_trimmed)
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


# ---- Lagged target -------------------------------------------------
y_lags_with_date <- create_lags(y_t%>% select(-sasdate), p_y, dates_yt)%>% select(-c(CPI_raw,CPI_t))
y_lags <- create_lags(y_t %>% select(-sasdate), p_y, dates_yt)%>% select(-c(sasdate, CPI_raw,CPI_t))

# ---- Lagged X (using stationary Xs) ------------------------------------------------------
X_t_lags <- create_lags(X_t, p_m, dates_Xt)

# ---- MARX & MAF ----
marx_data <- create_marx(X_t, max_lag = 12, dates_Xt)
maf_data <- create_maf(X_t, n_lags = 12, n_pcs = 2, dates_Xt)



# ---------------------- 3. Build Z-matrices (based on paper, every Z will include y-lags regardless + sasdate) ---------------


#F: Using factors only
#Z_F_stationary <- align_by_date(F_lags_stationary, y_lags)
#Z_F_raw <- align_by_date(F_lags_raw, y_lags)
Z_F_naked <- y_lags_with_date # for each naked, need to make 2 on ur end, one for Fraw and one for F stationary

#F-Level: Using factors + Levels
#Z_F_Level_stationary <- align_by_date(F_lags_stationary,y_t, y_lags)
#Z_F_Level_raw <- align_by_date(F_lags_raw, y_t, y_lags)
Z_F_Level_naked <- align_by_date(y_t, y_lags)

#X: Using stationary lagged Xs only
Z_X <- align_by_date(X_t_lags, y_lags)

#Using raw Xs only (not in paper)
Z_Ht <- align_by_date(X_t_raw_withdate, y_lags)

#X-MARX: Using stationary lagged Xs + MARX
Z_X_MARX <- align_by_date(X_t_lags, marx_data, y_lags)

#F-X-MARX: Using factors + stationary lagged Xs + MARX (Coulombe’s top performer for tree methods)
#Z_Fstationary_X_MARX_ <- align_by_date(F_lags_stationary, X_t_lags, marx_data, y_lags)
#Z_Fraw_X_MARX_ <- align_by_date(F_lags_raw, X_t_lags, marx_data, y_lags)
Z_F_X_MARX_naked <- align_by_date(X_t_lags, marx_data, y_lags)

#F-X-MAF: Using factors + stationary lagged Xs + MAF
#Z_Fstationary_X_MAF <- align_by_date(F_lags_stationary, X_t_lags, maf_data, y_lags)
#Z_Fraw_X_MAF <- align_by_date(F_lags_raw, X_t_lags, maf_data, y_lags)
Z_F_X_MAF_naked <- align_by_date(X_t_lags, maf_data, y_lags)

#X-MAF: Using stationary lagged Xs + MAF  
Z_X_MAF <- align_by_date(X_t_lags, maf_data, y_lags)

#F-X-MARX-Level: Using factors + stationary lagged Xs + MARX + Levels (raw Xs and Y_t)  
#Z_Fstationary_X_MARX_Level <- align_by_date(F_lags_stationary, X_t_lags, marx_data, y_lags, y_t)
#Z_Fraw_X_MARX_Level <- align_by_date(F_lags_raw, X_t_lags, marx_data, y_lags, y_t)
Z_F_X_MARX_Level_naked <- align_by_date(X_t_lags, marx_data, y_lags, y_t)


##SAVE AS RDS FILES
saveRDS(
  list(
    #Z_F_stationary = Z_F_stationary,
    #Z_F_raw = Z_F_raw,
    Z_F_naked = Z_F_naked,
    #Z_F_Level_stationary = Z_F_Level_stationary,
    #Z_F_Level_raw = Z_F_Level_raw,
    Z_F_Level_naked = Z_F_Level_naked,
    Z_X = Z_X,
    Z_Ht = Z_Ht,
    Z_X_MARX = Z_X_MARX,
    #Z_Fstationary_X_MARX_ = Z_Fstationary_X_MARX_,
    #Z_Fraw_X_MARX_ = Z_Fraw_X_MARX_,
    Z_F_X_MARX_naked = Z_F_X_MARX_naked,
    #Z_Fstationary_X_MAF = Z_Fstationary_X_MAF,
    #Z_Fraw_X_MAF = Z_Fraw_X_MAF,
    Z_F_X_MAF_naked = Z_F_X_MAF_naked,
    Z_X_MAF = Z_X_MAF,
    #Z_Fstationary_X_MARX_Level = Z_Fstationary_X_MARX_Level,
    #Z_Fraw_X_MARX_Level = Z_Fraw_X_MARX_Level,
    Z_F_X_MARX_Level_naked = Z_F_X_MARX_Level_naked
  ),
  file = "all_Z_matrices.rds"
)

message("✅ Static feature engineering complete. Ready for rolling PCA.")




# Rolling-window PCA function
#rolling_pca_lags <- function(X_full, dates_full, train_end_idx, max_lags = 4, var_threshold = 0.8, min_train = 50) {
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




  
  # Helper: align by date
 #align_by_date <- function(...) {
    dfs <- list(...)
    min_rows <- min(sapply(dfs, nrow))
    dfs_trimmed <- lapply(dfs, function(x) x[(nrow(x)-min_rows+1):nrow(x), , drop=FALSE])
    Reduce(function(x,y) cbind(x, y[ , setdiff(names(y), "sasdate")]), dfs_trimmed)
  }
  


saveRDS(y_t, file = "y_t.rds")
saveRDS(X_t, file = "X_t.rds")
