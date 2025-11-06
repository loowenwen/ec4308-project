library(dplyr)
library(glmnet)
library(doParallel)
library(foreach)
library(ggplot2)

# --- Read Data ---
Z_objects <- readRDS("all_Z_matrices.rds")
Z_names <- names(Z_objects)

# --- Order by matrix size (smallest → largest) ---
Z_order <- names(sort(sapply(Z_objects[Z_names], ncol)))
Z_names <- Z_order

# --- Load or initialize results ---
ridge_results_path <- "ridge_forecasts_allZ_aligned_fixed.rds"
lambda_results_path <- "ridge_lambda_allZ_aligned_fixed.rds"

all_ridge_results <- if (file.exists(ridge_results_path)) readRDS(ridge_results_path) else list()
all_lambda_results <- if (file.exists(lambda_results_path)) readRDS(lambda_results_path) else list()

processed_Z <- names(all_ridge_results)
Z_to_process <- setdiff(Z_names, processed_Z)

cat("Already processed:", processed_Z, "\n")
cat("To process (ordered smallest → largest):", Z_to_process, "\n")

if (length(Z_to_process) == 0) stop("All Z matrices processed.")

# --- Helper metrics ---
MSE <- function(pred, truth) mean((pred - truth)^2, na.rm = TRUE)
RMSE <- function(pred, truth) sqrt(MSE(pred, truth))
MAE <- function(pred, truth) mean(abs(pred - truth), na.rm = TRUE)

# --- Clean matrix ---
clean_Z_matrix <- function(Z_mat, target_col = "CPI_t") {
  if (!target_col %in% colnames(Z_mat)) stop(paste("Missing target column:", target_col))
  Z_mat <- Z_mat[!is.na(Z_mat[[target_col]]), ]
  pred_cols <- setdiff(colnames(Z_mat), c("sasdate", target_col))
  Z_mat <- Z_mat[!apply(Z_mat[, pred_cols], 1, function(x) any(is.na(x))), ]
  Z_mat <- Z_mat[, !duplicated(colnames(Z_mat))]
  rownames(Z_mat) <- NULL
  return(Z_mat)
}

# --- Time-series CV for λ (leakage-safe) ---
time_series_lambda_cv <- function(X, y, h = 1, n_oos = 66, window_size = 120, lambda_grid) {
  y_h <- dplyr::lead(y, h)
  valid_idx <- !is.na(y_h)
  X <- as.matrix(X[valid_idx, ])
  y <- y_h[valid_idx]
  
  n_total <- nrow(X)
  n_insample <- n_total - n_oos
  
  X_in <- X[1:n_insample, ]
  y_in <- y[1:n_insample]
  
  origins <- seq(window_size, n_insample - 1)
  
  mse_per_lambda <- foreach(lam = lambda_grid, .combine = c, .packages = "glmnet") %dopar% {
    fold_errors <- numeric(length(origins))
    for (k in seq_along(origins)) {
      end <- origins[k]
      start <- end - window_size + 1
      X_train <- X_in[start:end, ]
      y_train <- y_in[start:end]
      X_val <- matrix(X_in[end + 1, ], nrow = 1)
      y_val <- y_in[end + 1]
      
      fit <- glmnet(X_train, y_train, alpha = 0, lambda = lam, standardize = TRUE)
      pred <- as.numeric(predict(fit, newx = X_val, s = lam))
      fold_errors[k] <- (pred - y_val)^2
    }
    mean(fold_errors, na.rm = TRUE)
  }
  
  idx_min <- which.min(mse_per_lambda)
  lambda_min <- lambda_grid[idx_min]
  lambda_1se <- lambda_min
  list(lambda_min = lambda_min, lambda_1se = lambda_1se)
}

# --- Rolling Ridge Forecast (aligned & interpretable) ---
ridge_rolling_h <- function(X, y, h, nprev, window_size, lambda, lambda_grid = NULL, ret_errors = TRUE) {
  y_h <- dplyr::lead(y, h)
  valid_idx <- !is.na(y_h)
  X <- as.matrix(X[valid_idx, ])
  y <- y_h[valid_idx]
  n <- nrow(X)
  
  preds <- rep(NA, nprev)
  actuals <- rep(NA, nprev)
  
  for (i in seq_len(nprev)) {
    end <- n - nprev + i - 1
    start <- max(1, end - window_size + 1)
    X_train <- X[start:end, ]
    y_train <- y[start:end]
    X_test <- matrix(X[end + 1, ], nrow = 1)
    y_true <- y[end + 1]
    
  
    if (!is.null(lambda_grid)) {
      cv_res <- time_series_lambda_cv(X_train, y_train, h = 1, n_oos = 10, window_size = min(window_size, nrow(X_train) - 10), lambda_grid)
       lambda <- cv_res$lambda_min
     }
    
    fit <- glmnet(X_train, y_train, alpha = 0, lambda = lambda, standardize = TRUE)
    preds[i] <- as.numeric(predict(fit, newx = X_test, s = lambda))
    actuals[i] <- y_true
  }
  
  out <- list(
    preds = preds,
    real = actuals,
    rmse = RMSE(preds, actuals),
    mae = MAE(preds, actuals)
  )
  
  if (ret_errors) out$errors <- preds - actuals
  return(out)
}

# --- Parallel setup ---
n_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Running on", n_cores, "cores\n")

# --- Parameters ---
horizons <- c(1, 3, 6, 12)
n_oos <- 66
window_size <- 120
lambda_grid <- 10^seq(2, -3, length.out = 30)

# --- Main loop ---
for (z_name in Z_to_process) {
  cat("\n=== Ridge Forecast for", z_name, "===\n")
  
  tryCatch({
    Z_mat <- Z_objects[[z_name]] %>% arrange(sasdate)
    Z_mat <- clean_Z_matrix(Z_mat, "CPI_t")
    
    y_full <- Z_mat$CPI_t
    X_full <- Z_mat %>% select(-sasdate, -CPI_t)
    X_full <- model.matrix(~ ., data = X_full)[, -1]
    
    ridge_out <- list()
    lambda_out <- list()
    
    for (h in horizons) {
      cat("  Horizon:", h, "\n")
      
      # λ cross-validation (no leakage)
      cv_res <- time_series_lambda_cv(X_full, y_full, h, n_oos, window_size, lambda_grid)
      best_lambda <- cv_res$lambda_1se
      cat("    λ chosen:", signif(best_lambda, 3), "\n")
      lambda_out[[paste0("h", h)]] <- cv_res
      
      # Rolling forecast
      res <- ridge_rolling_h(X_full, y_full, h, nprev = n_oos, window_size = window_size, lambda = best_lambda)
      ridge_out[[paste0("h", h)]] <- res
    }
    
    all_ridge_results[[z_name]] <- ridge_out
    all_lambda_results[[z_name]] <- lambda_out
    
    saveRDS(all_ridge_results, ridge_results_path)
    saveRDS(all_lambda_results, lambda_results_path)
    cat("Saved results for:", z_name, "\n")
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

stopCluster(cl)
cat("\n✅ All Ridge forecasts completed (ordered smallest → largest Z matrices).\n")

# --- Reload results for plotting ---
ridge_results <- readRDS(ridge_results_path)
Z_objects <- readRDS("all_Z_matrices.rds")

# --- Plotting ---
out_dir <- "ridge_plots_clean_fixed"
if (!dir.exists(out_dir)) dir.create(out_dir)
horizons <- c("h1", "h3", "h6", "h12")

theme_white <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = "top",
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

for (z_name in names(ridge_results)) {
  cat("Plotting for:", z_name, "\n")
  
  ridge_out <- ridge_results[[z_name]]
  Z_mat <- Z_objects[[z_name]] %>% arrange(sasdate)
  
  for (h in horizons) {
    if (!h %in% names(ridge_out)) next
    
    preds <- ridge_out[[h]]$preds
    real <- ridge_out[[h]]$real
    nprev <- length(preds)
    sasdates <- tail(Z_mat$sasdate, nprev)
    
    df_plot <- data.frame(
      Date = sasdates,
      Actual = real,
      Predicted = preds
    )
    
    rmse_val <- round(ridge_out[[h]]$rmse, 4)
    mae_val  <- round(ridge_out[[h]]$mae, 4)
    
    p <- ggplot(df_plot, aes(x = Date)) +
      geom_line(aes(y = Actual, color = "Actual"), linewidth = 0.9) +
      geom_line(aes(y = Predicted, color = "Predicted"), linetype = "dashed", linewidth = 0.9) +
      scale_color_manual(values = c("Actual" = "#1A1A1A", "Predicted" = "#E63946")) +
      labs(
        title = paste("Actual vs Predicted —", z_name),
        subtitle = paste("H =", gsub("h", "", h), "| RMSE =", rmse_val, "| MAE =", mae_val),
        y = "CPI_t",
        x = "Date"
      ) +
      theme_white
    
    file_name <- file.path(out_dir, paste0(z_name, "_", h, "_plot.png"))
    ggsave(file_name, plot = p, width = 7, height = 4, dpi = 300)
  }
}

cat("\n✅ Clean plots with aligned forecasts saved to 'ridge_plots_clean_fixed' folder.\n")
