library(stats)     
library(zoo)   
library(tidyverse)
library(purrr)

# ---------------------- 1. Helper Functions -------------------------

#function to create lag variables
create_lags <- function(df, p, dates) {
  stopifnot(nrow(df) == length(dates))
  if (p == 0) return(tibble(sasdate = dates, df))
  
  lagged <- lapply(1:p, function(l) lag(df, l))
  lagged_df <- do.call(cbind, lagged)
  colnames(lagged_df) <- paste0(rep(colnames(df), each = p), "_L", 1:p)
  
  valid_rows <- (p + 1):nrow(df)
  tibble(sasdate = dates[valid_rows],
         as.data.frame(lagged_df[valid_rows, , drop = FALSE]))
}

#align by common date
align_by_date <- function(...) {
  dfs <- list(...)
  
  # 1. Find common dates across all input tibbles
  common_dates <- Reduce(intersect, lapply(dfs, `[[`, "sasdate"))
  
  # 2. Align and drop NA rows
  dfs_aligned <- lapply(dfs, function(x) {
    x %>%
      filter(sasdate %in% common_dates) %>%
      drop_na()
  })
  
  # 3. Merge everything into a single tibble by "sasdate"
  df_merged <- Reduce(function(x, y) left_join(x, y, by = "sasdate"), dfs_aligned)
  
  # 4. Return merged dataset
  return(df_merged)
}

#MARX Function (general idea: each column is the last p lags of each variable)
create_marx <- function(df, max_lag, dates) {          
  stopifnot(nrow(df) == length(dates))
  n <- nrow(df); k <- ncol(df)
  marx_matrix <- matrix(NA, nrow = n - max_lag, ncol = k * max_lag)
  
  for (p in 1:max_lag) {
    for (k_idx in 1:k) {
      col_idx <- (p - 1) * k + k_idx
      x <- as.matrix(df[, k_idx])
      # embed(x, p) gives p lags (no current value)
      x_embed <- embed(x, p)
      x_ma    <- rowMeans(x_embed, na.rm = TRUE) / p
      marx_matrix[, col_idx] <- tail(x_ma, n - max_lag)
    }
  }
  
  marx_df <- as.data.frame(marx_matrix)
  colnames(marx_df) <- paste0(rep(colnames(df), each = max_lag), "_MA", 1:max_lag)
  
  # return tibble with sasdate
  tibble(sasdate = tail(dates, n - max_lag), marx_df)
}

#MAF Function
create_maf <- function(df, n_lags = 12, n_pcs = 1, dates) {   # added dates
  stopifnot(nrow(df) == length(dates))
  n <- nrow(df); k <- ncol(df)
  maf_matrix <- matrix(NA, nrow = n, ncol = k * n_pcs)
  colnames(maf_matrix) <- paste0(rep(colnames(df), each = n_pcs), "_PC", 1:n_pcs)
  
  for (k_idx in 1:k){
    x <- as.matrix(df[, k_idx])
    lagged_mat <- embed(x, n_lags)               # n_lags columns = lags only
    lagged_mat <- na.omit(lagged_mat)
    n_rows <- nrow(lagged_mat)
    
    pca_res <- prcomp(lagged_mat, scale. = TRUE)
    pcs <- pca_res$x[, 1:n_pcs, drop = FALSE]
    
    maf_matrix[(n_lags:n)[1:n_rows],
               ((k_idx-1)*n_pcs + 1):(k_idx*n_pcs)] <- pcs
  }
  
  # keep only rows where we have PCs (after n_lags)
  valid_rows <- (n_lags + 1):n
  tibble(sasdate = dates[valid_rows],
         as.data.frame(maf_matrix[valid_rows, , drop = FALSE]))
}



# ---------------------- 2. (FOR TRAINING) Preparing factor lags (F), MARX, MAF, X lag (X_lag) -------------------------

p_y <- 4  #Number of lags for target 
p_f <- 4  #Number of lags for factors 
p_m <- 4 #Number of lags for X

X_t <- X_train %>% select(-sasdate) #from data_cleaning.r (use this)
X_train_raw <- fred_clean %>% filter(sasdate < split_date) 
X_t_raw <- X_train_raw %>% select(-sasdate, -CPIAUCSL)  #use this
train_dates_Xt <- X_train$sasdate #when using stationary Xs
train_dates_Xtraw <- X_train_raw$sasdate #when using non stationary Xs
train_dates_yt <- y_train$sasdate #when using yt

# ---- Lagged F (on stationary Xs) ------------------------------------------------------

#PCA on stationary Xs
pca_stationaryX_result <- prcomp(X_t, scale. = TRUE)
plot(pca_stationaryX_result, type = "l", main = "Scree Plot of PCA Factors (On Stationary Xs)", npcs = 125)
summary(pca_stationaryX_result)
num_factors <- 33 #explains 80% of variance


F_t_stationary <- as.data.frame(pca_stationaryX_result$x[, 1:num_factors]) #F_t_stationary = 33 factors
colnames(F_t_stationary) <- paste0("F", 1:num_factors)
F_lags_stationary <- create_lags(F_t_stationary, p_f, train_dates_Xt) #F lag using stationary Xs


# ---- Lagged F (on raw Xs) ------------------------------------------------------

#PCA on stationary Xs
pca_rawX_result <- prcomp(X_t_raw, scale. = TRUE)
plot(pca_rawX_result, type = "l", main = "Scree Plot of PCA Factors (On Raw Xs)", npcs = 125)
summary(pca_rawX_result)
num_factors_raw <- 4 #explains 88% of variance

F_t_raw <- as.data.frame(pca_rawX_result$x[, 1:num_factors_raw]) #F_t_raw = 4 factors
colnames(F_t_raw) <- paste0("F", 1:num_factors_raw)
F_lags_raw <- create_lags(F_t_raw, p_f, train_dates_Xtraw) #F lag using stationary Xs



# ---- y_t Target (CPI Level) ------------------------
y_t <- y_train

# ---- Lagged target -------------------------------------------------
y_lags <- create_lags(y_t %>% select(-sasdate), p_y, train_dates_yt)

# ---- Lagged X (using stationary Xs) ------------------------------------------------------
X_t_lags <- create_lags(X_t, p_m, train_dates_Xt)


# ---- MARX ----------------------------------------------------------
marx_data <- create_marx(X_t, max_lag = 12, train_dates_Xt)

# ---- MAF -----------------------------------------------------------
p_maf <- 12
n_pcs <- 2
maf_data  <- create_maf(X_t, n_lags = p_maf, n_pcs = n_pcs, train_dates_Xt)



# ---------------------- 3. (FOR TRAINING) Build Z-matrices (based on paper, every Z will include y-lags regardless + sasdate) ---------------


#F: Using factors only
Z_F_stationary <- align_by_date(F_lags_stationary, y_lags)
Z_F_raw <- align_by_date(F_lags_raw, y_lags)

#F-Level: Using factors + Levels
Z_F_stationary <- align_by_date(F_lags_stationary, y_t, y_lags)
Z_F_raw <- align_by_date(F_lags_raw, y_t, y_lags)

#X: Using stationary lagged Xs only
Z_X <- align_by_date(X_t_lags, y_lags)

#Using raw Xs only (not in paper)
Z_Ht <- align_by_date(X_train_raw, y_lags)

#X-MARX: Using stationary lagged Xs + MARX
Z_X_MARX <- align_by_date(X_t_lags, marx_data, y_lags)

#F-X-MARX: Using factors + stationary lagged Xs + MARX (Coulombe’s top performer for tree methods)
Z_Fstationary_X_MARX_ <- align_by_date(F_lags_stationary, X_t_lags, marx_data, y_lags)
Z_Fraw_X_MARX_ <- align_by_date(F_lags_raw, X_t_lags, marx_data, y_lags)

#F-X-MAF: Using factors + stationary lagged Xs + MAF
Z_Fstationary_X_MAF <- align_by_date(F_lags_stationary, X_t_lags, maf_data, y_lags)
Z_Fraw_X_MAF <- align_by_date(F_lags_raw, X_t_lags, maf_data, y_lags)

#X-MAF: Using stationary lagged Xs + MAF  
Z_X_MAF <- align_by_date(X_t_lags, maf_data, y_lags)

#F-X-MARX-Level: Using factors + stationary lagged Xs + MARX + Levels (raw Xs and Y_t)  
Z_Fstationary_X_MARX_Level <- align_by_date(F_lags_stationary, X_t_lags, maf_data, y_lags, y_t)
Z_Fraw_X_MARX_Level <- align_by_date(F_lags_raw, X_t_lags, maf_data, y_lags, y_t)




######################################CODES BELOW ARE FOR TEST SET##########################################################




# ---------------------- 4. (FOR TEST) Preparing factor lags (F), MARX, MAF, X lag (X_lag) -------------------------
X_test_raw <- fred_clean %>% filter(sasdate >= split_date) 
X_t_test <- X_test %>% select(-sasdate)  # from data_cleaning.r (stationary test set)
test_dates_Xt <- X_test$sasdate
test_dates_Xtraw <- X_test_raw$sasdate
test_dates_yt <- y_test$sasdate

# ---- Lagged F (on stationary Xs) ------------------------------------------------------
pca_stationaryX_result_test <- prcomp(X_t_test, scale. = TRUE)  # X_test_pca = test stationary predictors
plot(pca_stationaryX_result_test, type = "l", main = "Scree Plot (Test - Stationary Xs)", npcs = 125)
summary(pca_stationaryX_result_test)
num_factors_test <- 14
F_t_stationary_test <- as.data.frame(pca_stationaryX_result_test$x[, 1:num_factors_test])
colnames(F_t_stationary_test) <- paste0("F", 1:num_factors_test)
F_lags_stationary_test <- create_lags(F_t_stationary_test, p_f, test_dates_Xt)

# ---- Lagged F (on raw Xs) ------------------------------------------------------
pca_rawX_result_test <- prcomp((X_test_raw %>% select(-sasdate)), scale. = TRUE)  # X_test_raw = test raw predictors
plot(pca_rawX_result_test, type = "l", main = "Scree Plot (Test - Raw Xs)", npcs = 125)
summary(pca_rawX_result_test)
num_factors_test_raw <- 3
F_t_raw_test <- as.data.frame(pca_rawX_result_test$x[, 1:num_factors_test_raw])
colnames(F_t_raw_test) <- paste0("F", 1:num_factors_test_raw)
F_lags_raw_test <- create_lags(F_t_raw_test, p_f, test_dates_Xtraw)

# ---- y_t Target (CPI Level) ------------------------
y_t_test <- y_test

# ---- Lagged target -------------------------------------------------
y_lags_test <- create_lags(y_t_test %>% select(-sasdate), p_y, test_dates_yt)

# ---- Lagged X (using stationary Xs) ------------------------------------------------------
X_t_lags_test <- create_lags(X_t_test, p_m, test_dates_Xt)

# ---- MARX ----------------------------------------------------------
marx_data_test <- create_marx(X_t_test, max_lag = 12, test_dates_Xt)

# ---- MAF -----------------------------------------------------------
p_maf <- 12
n_pcs <- 2
maf_data_test <- create_maf(X_t_test, n_lags = p_maf, n_pcs = n_pcs, test_dates_Xt)

# ---------------------- 5. (FOR TEST) Build Z-matrices (same as train) ---------------
# F: Using factors only
Z_F_stationary_test <- align_by_date(F_lags_stationary_test, y_lags_test)
Z_F_raw_test <- align_by_date(F_lags_raw_test, y_lags_test)

# F-Level: Using factors + Levels
Z_F_stationary_test <- align_by_date(F_lags_stationary_test, y_t_test, y_lags_test)
Z_F_raw_test <- align_by_date(F_lags_raw_test, y_t_test, y_lags_test)

# X: Using stationary lagged Xs only
Z_X_test <- align_by_date(X_t_lags_test, y_lags_test)

# Using raw Xs only (not in paper)
Z_Ht_test <- align_by_date(X_test_raw, y_lags_test)

# X-MARX: Using stationary lagged Xs + MARX
Z_X_MARX_test <- align_by_date(X_t_lags_test, marx_data_test, y_lags_test)

# F-X-MARX: Using factors + stationary lagged Xs + MARX
Z_Fstationary_X_MARX_test <- align_by_date(F_lags_stationary_test, X_t_lags_test, marx_data_test, y_lags_test)
Z_Fraw_X_MARX_test <- align_by_date(F_lags_raw_test, X_t_lags_test, marx_data_test, y_lags_test)

# F-X-MAF: Using factors + stationary lagged Xs + MAF
Z_Fstationary_X_MAF_test <- align_by_date(F_lags_stationary_test, X_t_lags_test, maf_data_test, y_lags_test)
Z_Fraw_X_MAF_test <- align_by_date(F_lags_raw_test, X_t_lags_test, maf_data_test, y_lags_test)

# X-MAF: Using stationary lagged Xs + MAF
Z_X_MAF_test <- align_by_date(X_t_lags_test, maf_data_test, y_lags_test)

# F-X-MARX-Level: Using factors + stationary lagged Xs + MARX + Levels
Z_Fstationary_X_MARX_Level_test <- align_by_date(F_lags_stationary_test, X_t_lags_test, marx_data_test, y_lags_test, y_t_test)
Z_Fraw_X_MARX_Level_test <- align_by_date(F_lags_raw_test, X_t_lags_test, marx_data_test, y_lags_test, y_t_test)




##SAVE AS RDS FILES
saveRDS(
  list(
    # ---- TEST SET ----
    Z_F_stationary_test = Z_F_stationary_test,
    Z_F_raw_test = Z_F_raw_test,
    Z_X_test = Z_X_test,
    Z_Ht_test = Z_Ht_test,
    Z_X_MARX_test = Z_X_MARX_test,
    Z_Fstationary_X_MARX_test = Z_Fstationary_X_MARX_test,
    Z_Fraw_X_MARX_test = Z_Fraw_X_MARX_test,
    Z_Fstationary_X_MAF_test = Z_Fstationary_X_MAF_test,
    Z_Fraw_X_MAF_test = Z_Fraw_X_MAF_test,
    Z_X_MAF_test = Z_X_MAF_test,
    Z_Fstationary_X_MARX_Level_test = Z_Fstationary_X_MARX_Level_test,
    Z_Fraw_X_MARX_Level_test = Z_Fraw_X_MARX_Level_test,
    
    # ---- TRAIN SET ----
    Z_F_stationary = Z_F_stationary,
    Z_F_raw = Z_F_raw,
    Z_X = Z_X,
    Z_Ht = Z_Ht,
    Z_X_MARX = Z_X_MARX,
    Z_Fstationary_X_MARX_ = Z_Fstationary_X_MARX_,
    Z_Fraw_X_MARX_ = Z_Fraw_X_MARX_,
    Z_Fstationary_X_MAF = Z_Fstationary_X_MAF,
    Z_Fraw_X_MAF = Z_Fraw_X_MAF,
    Z_X_MAF = Z_X_MAF,
    Z_Fstationary_X_MARX_Level = Z_Fstationary_X_MARX_Level,
    Z_Fraw_X_MARX_Level = Z_Fraw_X_MARX_Level
  ),
  file = "all_Z_matrices.rds"
)

## HOW TO RUN: Z_objects <- readRDS("all_Z_matrices.rds") (list of all)







