
library(glmnet)
library(dplyr)
library(irlba)
library(foreach)
library(doParallel)


# ---------------- 1. PCA helper function ----------------
add_pca_factors <- function(X_train, X_test, n_pcs = 32, n_lags = 1) {
  pcs_basis <- prcomp(scale(X_train), center = TRUE, scale. = TRUE, rank. = n_pcs)
  
  pcs_train <- predict(pcs_basis, newdata = scale(X_train, center = pcs_basis$center, scale = pcs_basis$scale))
  pcs_train_df <- as.data.frame(pcs_train)
  colnames(pcs_train_df) <- paste0("PC", 1:n_pcs)
  
  pcs_test <- predict(pcs_basis, newdata = scale(X_test, center = pcs_basis$center, scale = pcs_basis$scale))
  pcs_test_df <- as.data.frame(pcs_test)
  colnames(pcs_test_df) <- paste0("PC", 1:n_pcs)
  
  # Add lagged PCs for training
  if (n_lags > 0) {
    for (lag_i in 1:n_lags) {
      lagged_train <- dplyr::lag(pcs_train_df, lag_i)
      colnames(lagged_train) <- paste0("PC", 1:n_pcs, "_lag", lag_i)
      pcs_train_df <- cbind(pcs_train_df, lagged_train)
    }
  }
  
  # Add lagged PCs for test
  if (n_lags > 0) {
    for (lag_i in 1:n_lags) {
      pcs_test_df[paste0("PC", 1:n_pcs, "_lag", lag_i)] <- 
        tail(pcs_train_df[, paste0("PC", 1:n_pcs)], lag_i)[lag_i, , drop = FALSE]
    }
  }
  
  list(train_pcs = pcs_train_df, test_pcs = pcs_test_df)
}

# ---------------- 2. Parameters ----------------
horizons <- c(1,3,6,12)
nprev <- 100
n_pcs <- 32
n_lags <- 4

save_path <- "elasticnet_results_100.rds"

# ---------------- 3. Load or initialize results ----------------
if (file.exists(save_path)) {
  results <- readRDS(save_path)
  cat("Loaded existing results from", save_path, "\n")
} else {
  results <- list()
  cat("Starting new results list\n")
}

# ---------------- 4. Load feature matrices ----------------
all_Z_list <- readRDS("all_Z_matrices.rds")
Z_list <- all_Z_list[!grepl("MAF", names(all_Z_list))]
Z_order <- names(sort(sapply(Z_list, function(z) ncol(z))))
Z_list <- Z_list[Z_order]

cat("Feature sets to run:\n")
print(names(Z_list))

# ---------------- 5. Setup parallel backend ----------------
n_cores <- parallel::detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Running Elastic Net using", n_cores, "cores\n")

# ---------------- 6. Common grids ----------------
alpha_grid <- seq(0, 1, length = 6)     # 0, 0.2, 0.4, 0.6, 0.8, 1
lambda_grid <- 10^seq(5, -2, length = 50)

# ======================================================
# NON-FACTOR (non-F) FEATURE SETS
# ======================================================
nonF_Z_list <- Z_list[!grepl("F", names(Z_list))]

for (z_name in names(nonF_Z_list)) {
  Z <- nonF_Z_list[[z_name]] %>% select(-sasdate)
  cat("\n=============================\nRunning for:", z_name, "\n")
  
  if (is.null(results[[z_name]])) results[[z_name]] <- list()
  
  # Align y_t and X_t
  if (nrow(y_t) > nrow(Z)) {
    y_aligned <- tail(y_t, nrow(Z))
    X_aligned <- tail(X_t, nrow(Z))
  } else {
    y_aligned <- y_t
    X_aligned <- X_t
  }
  
  for (h in horizons) {
    hname <- paste0("h", h)
    if (!is.null(results[[z_name]][[hname]])) {
      cat("  Skipping horizon", h, "(already done)\n")
      next
    }
    
    cat("  Horizon:", h, "\n")
    y_h <- dplyr::lead(y_aligned$CPI_t, h)
    valid_idx <- 1:(length(y_h) - h)
    
    preds <- foreach(i = seq_len(nprev), .combine = c,
                     .packages = c("glmnet", "dplyr")) %dopar% {
                       
                       train_start <- i
                       train_end <- length(y_h[valid_idx]) - nprev + i - 1
                       
                       Z_train <- Z[train_start:train_end, , drop = FALSE]
                       y_train <- y_h[train_start:train_end]
                       Z_test  <- Z[train_end + 1, , drop = FALSE]
                       
                       Z_train <- Z_train[, colSums(is.na(Z_train)) < nrow(Z_train), drop = FALSE]
                       Z_test  <- Z_test[, colnames(Z_train), drop = FALSE]
                       
                       valid_rows <- complete.cases(Z_train, y_train)
                       Z_train <- Z_train[valid_rows, , drop = FALSE]
                       y_train <- y_train[valid_rows]
                       
                       # --- Elastic Net tuning ---
                       best_mse <- Inf
                       best_alpha <- NA
                       best_model <- NULL
                       
                       for (a in alpha_grid) {
                         model_cv <- cv.glmnet(
                           as.matrix(Z_train), y_train,
                           alpha = a,
                           lambda = lambda_grid,
                           standardize = TRUE
                         )
                         if (min(model_cv$cvm) < best_mse) {
                           best_mse <- min(model_cv$cvm)
                           best_alpha <- a
                           best_model <- model_cv
                         }
                       }
                       
                       as.numeric(predict(best_model, newx = as.matrix(Z_test), s = "lambda.min"))
                     }
    
    real <- tail(y_h[valid_idx], nprev)
    rmse <- sqrt(mean((real - preds)^2))
    mae  <- mean(abs(real - preds))
    
    results[[z_name]][[hname]] <- list(pred = preds, errors = c(rmse = rmse, mae = mae))
    cat("    Horizon", h, "done. RMSE:", round(rmse, 4), "MAE:", round(mae, 4), "\n")
    saveRDS(results, save_path)
  }
  
  cat("*** Finished all horizons for", z_name, "***\n")
}

# ======================================================
# FACTOR-AUGMENTED (F) FEATURE SETS
# ======================================================
F_Z_list <- Z_list[grepl("F", names(Z_list))]
F_Z_list <- F_Z_list[4]


for (z_name in names(F_Z_list)) {
  Z <- F_Z_list[[z_name]] %>% select(-sasdate)
  cat("\n=============================\nRunning for:", z_name, "\n")
  
  if (is.null(results[[z_name]])) results[[z_name]] <- list()
  
  if (nrow(y_t) > nrow(Z)) {
    y_aligned <- tail(y_t, nrow(Z))
    X_aligned <- tail(X_t, nrow(Z))
  } else {
    y_aligned <- y_t
    X_aligned <- X_t
  }
  
  for (h in horizons) {
    hname <- paste0("h", h)
    if (!is.null(results[[z_name]][[hname]])) {
      cat("  Skipping horizon", h, "(already done)\n")
      next
    }
    
    cat("  Horizon:", h, "\n")
    y_h <- dplyr::lead(y_aligned$CPI_t, h)
    valid_idx <- 1:(length(y_h) - h)
    
    preds <- foreach(i = seq_len(nprev), .combine = c,
                     .packages = c("glmnet", "dplyr", "irlba")) %dopar% {
                       
                       train_start <- i
                       train_end <- length(y_h[valid_idx]) - nprev + i - 1
                       
                       Z_train_raw <- Z[train_start:train_end, ]
                       y_train <- y_h[train_start:train_end]
                       Z_test_raw <- Z[train_end + 1, , drop = FALSE]
                       
                       X_train_raw <- X_aligned[train_start:train_end, ]
                       X_test_raw  <- X_aligned[train_end + 1, , drop = FALSE]
                       
                       pcs_data <- add_pca_factors(X_train_raw, X_test_raw, n_pcs = n_pcs, n_lags = n_lags)
                       Z_train <- cbind(Z_train_raw, pcs_data$train_pcs)
                       Z_test  <- cbind(Z_test_raw, pcs_data$test_pcs)
                       
                       common_cols <- intersect(colnames(Z_train), colnames(Z_test))
                       Z_train <- Z_train[, common_cols, drop = FALSE]
                       Z_test  <- Z_test[, common_cols, drop = FALSE]
                       
                       Z_train <- Z_train[, colSums(is.na(Z_train)) < nrow(Z_train), drop = FALSE]
                       Z_test  <- Z_test[, colnames(Z_train), drop = FALSE]
                       
                       valid_rows <- complete.cases(Z_train, y_train)
                       Z_train <- Z_train[valid_rows, , drop = FALSE]
                       y_train <- y_train[valid_rows]
                       
                       # --- Elastic Net tuning ---
                       best_mse <- Inf
                       best_alpha <- NA
                       best_model <- NULL
                       
                       for (a in alpha_grid) {
                         model_cv <- cv.glmnet(
                           as.matrix(Z_train), y_train,
                           alpha = a,
                           lambda = lambda_grid,
                           standardize = TRUE
                         )
                         if (min(model_cv$cvm) < best_mse) {
                           best_mse <- min(model_cv$cvm)
                           best_alpha <- a
                           best_model <- model_cv
                         }
                       }
                       
                       as.numeric(predict(best_model, newx = as.matrix(Z_test), s = "lambda.min"))
                     }
    
    real <- tail(y_h[valid_idx], nprev)
    rmse <- sqrt(mean((real - preds)^2))
    mae  <- mean(abs(real - preds))
    
    results[[z_name]][[hname]] <- list(pred = preds, errors = c(rmse = rmse, mae = mae))
    cat("    Horizon", h, "done. RMSE:", round(rmse, 4), "MAE:", round(mae, 4), "\n")
    saveRDS(results, save_path)
  }
  
  cat("*** Finished all horizons for", z_name, "***\n")
}

# --- Clean up ---
stopCluster(cl)
cat("\nAll Elastic Net forecasts complete. Results saved to", save_path, "\n")
