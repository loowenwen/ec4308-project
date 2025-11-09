#ridge_results_100_sep <- lapply(ridge_results_100, function(model) {
  #lapply(model, function(horizon) {  list(
    #  pred = horizon$pred,
     # rmse = horizon$errors["rmse"],
     # mae  = horizon$errors["mae"]
#    )  }) })
#saveRDS(ridge_results_100, "ridge_results_100.rds")




# Purpose: Run ridge regression forecasts for both factor-augmented (F) and non-factor (non-F) datasets safely

library(glmnet)
library(dplyr)
library(irlba)
library(foreach)
library(doParallel)


add_pca_factors <- function(X_train, X_test, n_pcs = 32, n_lags = 1) {
  # Fit PCA on training
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
      # Take the last `lag_i` rows of the TRAINING set (before the forecast)
      # This ensures the test lag reflects the latest known values
      pcs_test_df[paste0("PC", 1:n_pcs, "_lag", lag_i)] <- 
        tail(pcs_train_df[, paste0("PC", 1:n_pcs)], lag_i)[lag_i, , drop = FALSE]
    }
  }
  
  list(train_pcs = pcs_train_df, test_pcs = pcs_test_df)
}

# --- Parameters ---
horizons <- c(1, 3, 6, 12)
nprev <- 100
n_pcs <- 32
n_lags <- 4

save_path <- "ridge_results_100.rds"

# ---------------- 2. Load or initialize results ----------------
if (file.exists(save_path)) {
  results <- readRDS(save_path)
  cat("Loaded existing results from", save_path, "\n")
} else {
  results <- list()
  cat("Starting new results list\n")
}

# ---------------- 3. Load feature matrices ----------------
all_Z_list <- readRDS("all_Z_matrices.rds")

# Filter: remove any matrices with "MAF" in name
Z_list <- all_Z_list[!grepl("MAF", names(all_Z_list))]

# Keep everything else (F, X, MARX, etc.)
Z_order <- names(sort(sapply(Z_list, function(z) ncol(z))))
Z_list <- Z_list[Z_order]

cat("Feature sets to run:\n")
print(names(Z_list))

# ---------------- 4. Setup parallel backend ----------------
n_cores <- parallel::detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Running Ridge regression using", n_cores, "cores\n")



nonF_Z_list <- Z_list[!grepl("F", names(Z_list))]

for (z_name in names(nonF_Z_list)) {
  Z <- nonF_Z_list[[z_name]] %>% select(-sasdate)
  cat("\n=============================\nRunning for:", z_name, "\n")
  
  if (is.null(results[[z_name]])) results[[z_name]] <- list()
  
  # Align y_t and X_t to Z’s row count
  if (nrow(y_t) > nrow(Z)) {
    y_aligned <- tail(y_t, nrow(Z))
    X_aligned <- tail(X_t, nrow(Z))  # may not be needed for non-F
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
                       
                       # Drop any fully NA columns
                       Z_train <- Z_train[, colSums(is.na(Z_train)) < nrow(Z_train), drop = FALSE]
                       Z_test  <- Z_test[, colnames(Z_train), drop = FALSE]
                       
                       # Drop rows with NA in training set
                       valid_rows <- complete.cases(Z_train, y_train)
                       Z_train <- Z_train[valid_rows, , drop = FALSE]
                       y_train <- y_train[valid_rows]
                       
                       # Ridge regression
                       lambda_grid <- 10^seq(5, -2, length = 50)
                       model <- cv.glmnet(as.matrix(Z_train), y_train, alpha = 0, lambda = lambda_grid, standardize = TRUE)
                       as.numeric(predict(model, newx = as.matrix(Z_test), s = "lambda.min"))
                     }
    
    # Compute RMSE/MAE
    real <- tail(y_h[valid_idx], nprev)
    rmse <- sqrt(mean((real - preds)^2))
    mae  <- mean(abs(real - preds))
    
    results[[z_name]][[hname]] <- list(pred = preds, errors = c(rmse = rmse, mae = mae))
    cat("    Horizon", h, "done. RMSE:", round(rmse, 4), "MAE:", round(mae, 4), "\n")
    saveRDS(results, save_path)
  }
  
  cat("*** Finished all horizons for", z_name, "***\n")
}



###DOING THE Zs with Fs


F_Z_list <- Z_list[grepl("F", names(Z_list))]

for (z_name in names(F_Z_list)) {
  Z <- F_Z_list[[z_name]] %>% select(-sasdate)
  cat("\n=============================\nRunning for:", z_name, "\n")
  
  if (is.null(results[[z_name]])) results[[z_name]] <- list()
  
  # Align y_t and X_t to Z’s row count
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
                       
                       # === Compute PCA factors and add lagged PCs ===
                       pcs_data <- add_pca_factors(X_train_raw, X_test_raw, n_pcs = n_pcs, n_lags = n_lags)
                       Z_train <- cbind(Z_train_raw, pcs_data$train_pcs)
                       Z_test  <- cbind(Z_test_raw, pcs_data$test_pcs)
                       
                       # Align columns
                       common_cols <- intersect(colnames(Z_train), colnames(Z_test))
                       Z_train <- Z_train[, common_cols, drop = FALSE]
                       Z_test  <- Z_test[, common_cols, drop = FALSE]
                       
                       # Drop fully NA columns
                       Z_train <- Z_train[, colSums(is.na(Z_train)) < nrow(Z_train), drop = FALSE]
                       Z_test  <- Z_test[, colnames(Z_train), drop = FALSE]
                       
                       # Drop rows with NA
                       valid_rows <- complete.cases(Z_train, y_train)
                       Z_train <- Z_train[valid_rows, , drop = FALSE]
                       y_train <- y_train[valid_rows]
                       
                       # Ridge regression
                       lambda_grid <- 10^seq(5, -2, length = 50)
                       model <- cv.glmnet(as.matrix(Z_train), y_train, alpha = 0, lambda = lambda_grid, standardize = TRUE)
                       as.numeric(predict(model, newx = as.matrix(Z_test), s = "lambda.min"))
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



library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

# Load saved results (if not already)
results <- readRDS("ridge_results_100.rds")

# Assume y_t has full CPI series
actual_full <- y_t$CPI_t
dates_full  <- y_t$sasdate
nprev <- 100

# --- Helper: extract all forecasts for all horizons ---
extract_F_forecasts <- function(results, y_t, horizons, nprev) {
  plot_list <- list()  # initialize inside function
  
  # Filter only factor-augmented models (names containing "F")
  F_names <- names(results)[grepl("F", names(results))]
  
  for (z_name in F_names) {
    for (h in horizons) {
      hname <- paste0("h", h)
      
      preds <- results[[z_name]][[hname]]$pred
      
      # Lead y_t by h to align actuals with predictions
      actuals <- dplyr::lead(y_t$CPI_t, h)
      
      # Only take the portion corresponding to your nprev predictions
      actuals <- tail(actuals, nprev)
      
      plot_list[[paste(z_name, h, sep = "_")]] <- data.frame(
        time = tail(y_t$date, nprev),  # or the corresponding time index
        series = rep(c("actual", "pred"), each = length(preds)),
        value = c(actuals, preds),
        horizon = h,
        model = z_name
      )
    }
  }
  
  # Combine all into one dataframe
  plot_long <- dplyr::bind_rows(plot_list)
  return(plot_long)
}


ggplot(plot_long, aes(x = time, y = value, color = model, linetype = series, group = interaction(model, series))) +
  geom_line(alpha = 0.9, linewidth = 0.6) +
  facet_wrap(~ horizon, scales = "free_y", ncol = 2) +
  scale_color_brewer(palette = "Set1") +  # each Z gets a distinct color
  scale_linetype_manual(values = c("actual" = "solid", "pred" = "dashed")) +
  labs(
    title = "CPI Forecasts: Factor-Augmented Ridge Models",
    subtitle = "Actual vs Predicted across all horizons",
    x = "Time", y = "CPI",
    color = "Model", linetype = "Series"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )


plot_df <- extract_all_forecasts(results)

# --- Combine into tidy long format ---
plot_long <- plot_df %>%
  pivot_longer(cols = c(actual, pred), names_to = "series", values_to = "value")

# --- Plot all horizons and model types ---
ggplot(plot_long, aes(x = time, y = value, color = series, group = interaction(model, series))) +
  geom_line(alpha = 0.9, linewidth = 0.2) +
  facet_wrap(~ horizon, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("actual" = "black", "pred" = "#0072B2")) +
  labs(
    title = "CPI Forecasts: Factor-Augmented Ridge Models",
    subtitle = "Actual vs Predicted across all horizons",
    x = "Time", y = "CPI",
    color = "Series"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )



F_names <- names(ridge_results_100)[grepl("F", names(ridge_results_100))]

# Extract RMSEs
F_rmse_list <- lapply(F_names, function(z_name) {
  res_list <- ridge_results_100[[z_name]]
  sapply(res_list, function(h) h$errors["rmse"])
})

# Combine into a data frame
F_rmse_df <- do.call(rbind, lapply(seq_along(F_rmse_list), function(i) {
  data.frame(
    Z_name = F_names[i],
    Horizon = names(F_rmse_list[[i]]),
    RMSE = as.numeric(F_rmse_list[[i]])
  )
}))


##SANITY CHECK add_pca_factors###
#i <- 1
#train_start <- i
#train_end <- nrow(X_t) - nprev + i - 1

# Slice your X aligned to Z (or X_t)
#X_train_raw <- X_t[train_start:train_end, , drop = FALSE]
#X_test_raw  <- X_t[train_end + 1, , drop = FALSE]

# Now call the PCA function
#pcs_data <- add_pca_factors(X_train_raw, X_test_raw, n_pcs = 32, n_lags = 4)

# Check the results
#dim(pcs_data$train_pcs)
dim(pcs_data$test_pcs)
colnames(pcs_data$train_pcs)
colnames(pcs_data$test_pcs)