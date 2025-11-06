library(stats)
library(zoo)
library(tidyverse)
library(purrr)

# ---------------------- 1. Helper Functions -------------------------

# Create lag variables
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






