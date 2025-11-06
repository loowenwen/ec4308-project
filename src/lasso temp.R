library(hdm)
library(dplyr)
library(foreach)
library(doParallel)

horizons <- c(1,3,6,12)
save_path <- "rlasso_bic_results.rds"

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
nprev <- 66
Z_order <- names(sort(sapply(Z_list, function(z) ncol(z))))
Z_list <- Z_list[Z_order]

# Small lambda grid for BIC
lambda_grid <- c(0.01, 0.1, 1)

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
    out <- foreach(i = seq_len(nprev), .combine = 'c', .packages = c("hdm")) %dopar% {
      train_start <- i
      train_end <- length(y_h[valid_idx]) - nprev + i - 1
      Z_train <- Z[train_start:train_end, ]
      y_train <- y_h[train_start:train_end]
      Z_test <- Z[train_end + 1, , drop = FALSE]
      
      # BIC selection
      bic_vals <- sapply(lambda_grid, function(lambda) {
        fit <- rlasso(Z_train, y_train, post = FALSE, lambda = lambda)
        rss <- sum((y_train - predict(fit, newdata = Z_train))^2)
        df <- sum(coef(fit) != 0)
        n <- length(y_train)
        rss + log(n) * df
      })
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


# plot
