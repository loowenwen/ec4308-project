library(hdm)
library(dplyr)
library(foreach)
library(doParallel)
library(irlba)


horizons <- c(1,3,6,12)
save_path <- "rlasso_bic_results_100.rds"

# Load previous progress if available
if (file.exists(save_path)) {
  results <- readRDS(save_path)
  cat("Loaded existing results from", save_path, "\n")
} else {
  results <- list()
  cat("Starting new results list\n")
}

# Load and reorder Z matrices
all_Z_list <- readRDS("all_Z_matrices.rds")
Z_list <- all_Z_list
y_t <- y_t
nprev <- 100
Z_order <- names(sort(sapply(Z_list, function(z) ncol(z))))
Z_list <- Z_list[Z_order]

# Filter to keep only Zs without "F" in their names
Z_list <- Z_list[!grepl("F", names(Z_list))]

# Small lambda grid for BIC
#lambda_grid <- c(0.01, 0.1, 1)
lambda_grid <- 10^seq(-6,1,length.out=30)

# Parallel setup
n_cores <- parallel::detectCores() - 3
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores\n")

# Main loop
for (z_name in names(Z_list)) {
  Z_full <- Z_list[[z_name]]
  Z <- Z_full %>% select(-sasdate)
  cat("\n=============================\n")
  cat("Running for:", z_name, "\n")
  
  if (is.null(results[[z_name]])) results[[z_name]] <- list()
  
  y_aligned <- tail(y_t, nrow(Z))
  
  # Print first 3 dates for diagnostics
  cat("  First 3 Z dates:", paste(as.character(Z_full$sasdate[1:3]), collapse = ", "), "\n")
  cat("  First 3 y dates:", paste(as.character(y_aligned$sasdate[1:3]), collapse = ", "), "\n")
  
  for (h in horizons) {
    if (!is.null(results[[z_name]][[paste0("h", h)]])) {
      cat("  Skipping horizon", h, "(already done)\n")
      next
    }
    
    cat("  Horizon:", h, "\n")
    y_h <- dplyr::lead(y_aligned$CPI_t, h)
    valid_idx <- 1:(length(y_h) - h)
    
    # Parallel rolling window
    out <- foreach(i = seq_len(nprev), .combine = 'c', .packages = c("glmnet")) %dopar% {
      train_start <- i
      train_end <- length(y_h[valid_idx]) - nprev + i - 1
      Z_train <- as.matrix(Z[train_start:train_end, ])
      y_train <- y_h[train_start:train_end]
      Z_test <- as.matrix(Z[train_end + 1, , drop = FALSE])
      
      # Fit glmnet for the whole lambda grid at once
      fit <- glmnet(Z_train, y_train, alpha = 1, standardize = TRUE, intercept = TRUE, lambda = lambda_grid)
      
      # Compute BIC for each lambda
      y_hat <- predict(fit, Z_train)
      rss <- colSums((y_train - y_hat)^2)
      df <- fit$df  # non-zero coefficients per lambda
      n <- length(y_train)
      bic_vals <- rss + log(n) * df
      
      # Select lambda with minimum BIC
      best_lambda <- lambda_grid[which.min(bic_vals)]
      model <- rlasso(Z_train, y_train, post = FALSE, lambda = best_lambda)
      predict(model, newdata = Z_test)
    }
    
    results[[z_name]][[paste0("h", h)]] <- list(
      pred = out,
      errors = c(
        rmse = sqrt(mean((tail(y_h[valid_idx], nprev) - out)^2)),
        mae  = mean(abs(tail(y_h[valid_idx], nprev) - out))
      )
    )
    
    cat("    Horizon", h, "done.\n")
    saveRDS(results, save_path)
    cat("  Progress saved to", save_path, "\n")
  }
  
  cat("*** Finished all horizons for", z_name, "***\n\n")
}

stopCluster(cl)

# Summary
cat("\nSummary of results:\n")
for (z_name in names(results)) {
  cat("Results for:", z_name, "\n")
  for (h in horizons) {
    if (!is.null(results[[z_name]][[paste0("h", h)]])) {
      err <- results[[z_name]][[paste0("h", h)]]$errors
      cat("  Horizon", h, ": RMSE =", round(err["rmse"], 4), ", MAE =", round(err["mae"], 4), "\n")
    }
  }
}




# ====== adding PCs ===========================
library(dplyr)
library(irlba)

# Helper: augment Z with PCA factors and their lags
add_pca_factors <- function(Z, n_pcs = 32, n_lags = 1) {
  # Compute top PCs
  pcs <- get_pcs(Z, n_pcs = n_pcs)
  pcs_df <- as.data.frame(pcs)
  colnames(pcs_df) <- paste0("PC", 1:n_pcs)  # explicit names for the PCs themselves
  
  # Add lagged versions of PCs
  if (n_lags > 0) {
    for (lag_i in 1:n_lags) {
      lagged_pcs <- dplyr::lag(pcs_df, lag_i)
      colnames(lagged_pcs) <- paste0("PC", 1:n_pcs, "_lag", lag_i)
      pcs_df <- cbind(pcs_df, lagged_pcs)
    }
  }
  
  # Combine original Z and PCs (+ lags)
  Z_aug <- cbind(Z, pcs_df)
  
  return(Z_aug)
}

# --- Helper function: compute PCs and add lags ---
get_pcs <- function(Z, n_pcs = 32) {
  # Handles potential collinearity or rank-deficiency better than prcomp
  svd_res <- irlba(scale(Z), nv = n_pcs)
  return(svd_res$u %*% diag(svd_res$d))
}

add_pca_factors <- function(Z, n_pcs = 32, n_lags = 1) {
  pcs <- get_pcs(Z, n_pcs = n_pcs)
  pcs_df <- as.data.frame(pcs)
  colnames(pcs_df) <- paste0("PC", 1:n_pcs)
  
  # Add lagged PCs
  if (n_lags > 0) {
    for (lag_i in 1:n_lags) {
      lagged_pcs <- dplyr::lag(pcs_df, lag_i)
      colnames(lagged_pcs) <- paste0("PC", 1:n_pcs, "_lag", lag_i)
      pcs_df <- cbind(pcs_df, lagged_pcs)
    }
  }
  
  # Combine with original design matrix
  Z_aug <- cbind(Z, pcs_df)
  return(Z_aug)
}



save_path <- "rlasso_bic_results_100.rds"

# --- Main rolling window loop ---
horizons <- c(1, 3, 6, 12)
nprev <- 100
n_pcs <- 32
n_lags <- 4
save_path <- "rlasso_bic_results_100.rds"

if (file.exists(save_path)) {
  results <- readRDS(save_path)
  cat("Loaded existing results from", save_path, "\n")
} else {
  results <- list()
}



Z_list <- all_Z_list[grepl("F", names(all_Z_list))]  # only F-type matrices

for (z_name in names(Z_list)) {
  Z_full <- Z_list[[z_name]]
  Z <- Z_full %>% select(-sasdate)
  cat("\n=============================\nRunning for:", z_name, "\n")
  
  if (is.null(results[[z_name]])) results[[z_name]] <- list()
  
  # Align y
  if (nrow(y_t) > nrow(Z)) {
    y_aligned <- tail(y_t, nrow(Z))
  } else {
    y_aligned <- y_t
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
    
    # --- Parallel rolling window ---
    preds <- foreach(i = seq_len(nprev), .combine = c, .packages = c("hdm", "dplyr", "irlba")) %dopar% {
      train_start <- i
      train_end <- length(y_h[valid_idx]) - nprev + i - 1
      Z_train_raw <- Z[train_start:train_end, ]
      y_train <- y_h[train_start:train_end]
      Z_test_raw <- Z[train_end + 1, , drop = FALSE]
      
      # === Add PCA factors dynamically ===
      Z_train <- add_pca_factors(Z_train_raw, n_pcs = n_pcs, n_lags = n_lags)
      
      # Apply same PCA transformation to test (project test obs on training PCs)
      pcs_basis <- prcomp(scale(Z_train_raw), center = TRUE, scale. = TRUE, rank. = n_pcs)
      pcs_test <- predict(pcs_basis, newdata = scale(Z_test_raw, center = attr(pcs_basis$center, "scaled:center"), scale = attr(pcs_basis$scale, "scaled:scale")))
      pcs_test_df <- as.data.frame(pcs_test)
      colnames(pcs_test_df) <- paste0("PC", 1:n_pcs)
      
      # Add lags to test (fill NAs with last available value)
      if (n_lags > 0) {
        for (lag_i in 1:n_lags) {
          lagged_train <- tail(Z_train[, paste0("PC", 1:n_pcs)], lag_i)
          lagged_names <- paste0("PC", 1:n_pcs, "_lag", lag_i)
          pcs_test_df[lagged_names] <- lagged_train[1, ]
        }
      }
      
      # Combine with original Z_test
      Z_test <- cbind(Z_test_raw, pcs_test_df)
      
      # Plug-in LASSO
      model <- rlasso(Z_train, y_train, post = TRUE)
      predict(model, newdata = Z_test)
    }
    
    # --- Compute RMSE and MAE ---
    real <- tail(y_h[valid_idx], nprev)
    rmse <- sqrt(mean((real - preds)^2))
    mae <- mean(abs(real - preds))
    
    results[[z_name]][[hname]] <- list(pred = preds, errors = c(rmse = rmse, mae = mae))
    cat("    Horizon", h, "done. RMSE:", round(rmse, 4), "MAE:", round(mae, 4), "\n")
    saveRDS(results, save_path)
  }
  
  cat("*** Finished all horizons for", z_name, "***\n")
}

stopCluster(cl)


# plot
