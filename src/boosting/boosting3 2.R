library(gbm)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tibble)
library(parallel)

set.seed(123)

# Load Data 
y_train <- read.csv("../../data/y_train.csv")
y_test  <- read.csv("../../data/y_test.csv")
all_Z_matrices <- readRDS("~/Downloads/Y4:S1 NUS 2025:2026/EC4308/ec4308-project/src/all_Z_matrices.rds")

y_var <- bind_rows(y_train, y_test) %>%
  select(CPI_t, sasdate) %>%
  mutate(sasdate = ymd(sasdate),
         target_var = CPI_t) %>%
  select(-CPI_t)

y_train <- y_var %>%
  slice(1:(n() - 100))

y_test <- y_var %>%
  slice((n()-100)+1: n())

#--------------------------
# Helper Functions
#--------------------------
data_adjusted_time_horizon <- function(data_matrix, time_horizon){
  data_matrix %>%
    mutate(target_var = dplyr::lead(target_var, time_horizon),
           sasdate = dplyr::lead(sasdate, time_horizon)) %>%
    na.omit()
}

# Robust GBM training + validation (logs progress)
boosting_for_one_parameter <- function(train_boosting, val_boosting, n_trees, time_horizon, dataset_name) {
  if (nrow(val_boosting) == 0) return(list(mse = Inf, fitted_model = NULL))
  
  mse_vec <- numeric(nrow(val_boosting))
  current_train <- as_tibble(train_boosting)
  last_model <- NULL
  
  for (i in seq_len(nrow(val_boosting))) {
    # try to catch errors in fitting
    fit_res <- tryCatch({
      model <- gbm(target_var ~ ., data = current_train, distribution = "gaussian",
                   bag.fraction = 0.5, interaction.depth = 2, n.trees = n_trees,
                   shrinkage = 0.01, verbose = FALSE)
      list(success = TRUE, model = model)
    }, error = function(e) {
      list(success = FALSE, error = conditionMessage(e), model = NULL)
    })
    
    if (!fit_res$success) {
      cat(Sys.time(), " | ERROR fitting GBM:", fit_res$error,
          " | dataset:", dataset_name, " | horizon:", time_horizon, " | n_trees:", n_trees, "\n",
          file = paste0("progress_", dataset_name, "_h", time_horizon, ".log"), append = TRUE)
      # assign NA to the remaining mse entries and return
      mse_vec[i:length(mse_vec)] <- NA_real_
      return(list(mse = Inf, fitted_model = NULL))
    }
    
    last_model <- fit_res$model
    
    newdata_row <- val_boosting[i, setdiff(names(val_boosting), "target_var"), drop = FALSE] %>% as_tibble()
    pred <- predict(last_model, newdata = newdata_row, n.trees = n_trees)
    actual <- val_boosting$target_var[i]
    mse_vec[i] <- (pred - actual)^2
    
    # Log validation progress to file
    cat(Sys.time(),
        " | Val iter:", i,
        " | dataset:", dataset_name,
        " | horizon:", time_horizon,
        " | n_trees:", n_trees,
        " | pred:", signif(pred, 6),
        " | actual:", signif(actual, 6),
        "\n",
        file = paste0("progress_", dataset_name, "_h", time_horizon, ".log"),
        append = TRUE)
    
    # Safe rolling-origin update: drop oldest row, add observed val row
    if (nrow(current_train) > 1) {
      current_train <- bind_rows(tail(current_train, nrow(current_train) - 1) %>% as_tibble(), val_boosting[i, , drop = FALSE] %>% as_tibble())
    } else {
      current_train <- val_boosting[i, , drop = FALSE] %>% as_tibble()
    }
  }
  
  rmse <- sqrt(mean(mse_vec, na.rm = TRUE))
  list(mse = rmse, fitted_model = last_model)
}

# Finds best parameter (n_trees) by running boosting_for_one_parameter
gradient_boosting_best_model <- function(data_matrix_boost, parameters_boost, time_horizon, dataset_name) {
  train_boosting <- data_matrix_boost %>%
    head(nrow(data_matrix_boost) - 110) %>%
    select(-sasdate) %>%
    na.omit() %>%
    as_tibble()
  
  val_boosting <- data_matrix_boost %>%
    select(-sasdate) %>%
    tail(110) %>%
    as_tibble()
  
  best_mse <- Inf
  best_model <- NULL
  best_param <- NA
  
  for (n_trees in parameters_boost) {
    res <- boosting_for_one_parameter(train_boosting, val_boosting, n_trees, time_horizon, dataset_name)
    if (!is.null(res$fitted_model) && is.finite(res$mse) && res$mse < best_mse) {
      best_mse <- res$mse
      best_model <- res$fitted_model
      best_param <- n_trees
    }
  }
  
  list(best_model = best_model, best_param = best_param, best_mse = best_mse)
}

# Predict across rolling test set for one horizon; log predictions and save per-horizon file
predict_test_set_per_time_horizon <- function(data_matrix_boosting, parameters_boost, time_horizon, dataset_name, test_start_date = as.Date("2017-03-01")) {
  # prepare matrix and align the target for horizon
  data_matrix_boosting <- data_matrix_boosting %>%
    mutate(sasdate = ymd(sasdate)) %>%
    left_join(y_var, by = "sasdate") %>%
    data_adjusted_time_horizon(time_horizon)
  
  test_set <- data_matrix_boosting %>% filter(sasdate >= test_start_date) %>% as_tibble()
  train_set <- data_matrix_boosting %>% filter(sasdate < test_start_date) %>% as_tibble()
  
  if (nrow(test_set) == 0) {
    stop("No test rows found for given test_start_date")
  }
  
  # initialize preds tibble with fixed columns
  preds_df <- tibble(
    dataset = character(),
    time_horizon = numeric(),
    date_predicted = as.Date(character()),
    predicted_y = numeric(),
    actual_y = numeric(),
    rmse_val = numeric(), 
    best_param = numeric(),
  )
  
  initial_train_size <- nrow(train_set)
  
  # iterate through test rows and do rolling-origin updates
  for (k in seq_len(nrow(test_set))) {
    date_y_t <- test_set$sasdate[k]
    
    # get best model on current train_set
    best_model_info <- gradient_boosting_best_model(train_set, parameters_boost, time_horizon, dataset_name)
    model <- best_model_info$best_model
    n_trees <- best_model_info$best_param
    rmse_val <- best_model_info$best_mse
    
    # If no model returned, log and skip this k
    if (is.null(model)) {
      cat(Sys.time(), " | No model found for dataset:", dataset_name, "horizon:", time_horizon, "at date:", as.character(date_y_t), "\n",
          file = paste0("progress_", dataset_name, "_h", time_horizon, ".log"), append = TRUE)
      # Expand training set safely and continue
      if (nrow(train_set) > 1) {
        train_set <- bind_rows(tail(train_set, nrow(train_set) - 1) %>% as_tibble(), test_set[k, , drop = FALSE] %>% as_tibble())
      } else {
        train_set <- test_set[k, , drop = FALSE] %>% as_tibble()
      }
      next
    }
    
    # predict for current test row
    newdata_row <- test_set[k, setdiff(names(test_set), "target_var"), drop = FALSE] %>% as_tibble()
    pred_y <- predict(model, newdata = newdata_row, n.trees = n_trees)
    actual_y <- test_set$target_var[k]
    
    # Log the prediction
    cat(Sys.time(),
        " | Predicted on:", as.character(date_y_t),
        " | Horizon:", time_horizon,
        " | Dataset:", dataset_name,
        " | Pred:", signif(pred_y, 6),
        " | Actual:", signif(actual_y, 6),
        " | RMSE:", signif(rmse_val, 6),
        " | Parameter Choice:", signif(n_trees, 6),
        "\n",
        file = paste0("progress_", dataset_name, "_h", time_horizon, ".log"),
        append = TRUE)
    
    # add to preds_df
    preds_df <- add_row(preds_df,
                        dataset = dataset_name,
                        time_horizon = time_horizon,
                        date_predicted = date_y_t,
                        predicted_y = pred_y,
                        actual_y = actual_y,
                        rmse_val = rmse_val, 
                        best_param = n_trees)
    save_path <- paste0("predictions_", dataset_name, "_h", time_horizon, ".rds")
    saveRDS(preds_df, file = save_path)
    
    # Expand training set for rolling-origin safely
    if (nrow(train_set) > 1) {
      train_set <- bind_rows(tail(train_set, nrow(train_set) - 1) %>% as_tibble(), test_set[k, , drop = FALSE] %>% as_tibble())
    } else {
      train_set <- test_set[k, , drop = FALSE] %>% as_tibble()
    }
  }
  
  # Save per-horizon result immediately (caller may also save)
  
  cat(Sys.time(), " | Saved predictions to:", save_path, "\n", file = paste0("progress_", dataset_name, "_h", time_horizon, ".log"), append = TRUE)
  
  preds_df
}

####################
# Parameter Choice #
####################

## Best Z_X is found to be 200 - 300 for all horizons using OOB 
data_matrix_boosting <- all_Z_matrices[["Z_X"]] %>%
  mutate(sasdate = ymd(sasdate)) %>%
  left_join(y_var, by = "sasdate") %>%
  data_adjusted_time_horizon(12)
train <- data_matrix_boosting %>%
  na.omit() %>%
  slice(1:(n() - 100)) %>%
  select(-sasdate)
boostfit = gbm(target_var~.,data=train,distribution='gaussian',bag.fraction = .5,
               interaction.depth=2,n.trees=10000,shrinkage=.01)
best = gbm.perf(boostfit, method="OOB")


parameters_boost <- c(300, 500, 800)      
datasets<- c("Z_X")           
time_horizons <- c(1, 3, 6, 12)

#--------------------------
# Parallel Execution (each horizon runs in its own process)
#--------------------------
num_cores <- max(1, detectCores() - 1)
all_results <- list()

for (dataset_name in datasets) {
  data_matrix <- all_Z_matrices[[dataset_name]]
  
  # parallel across horizons
  results_list <- mclapply(time_horizons, function(h) {
    message("Starting dataset=", dataset_name, " horizon=", h)
    predict_test_set_per_time_horizon(data_matrix, parameters_boost, h, dataset_name)
  }, mc.cores = min(length(time_horizons), num_cores))
  
  # combine results for this dataset (each element is a tibble)
  all_results[[dataset_name]] <- bind_rows(results_list)
}

# Combine everything into one master tibble (optional)
final_preds <- bind_rows(all_results)

# Save combined results
saveRDS(final_preds, file = "predictions_parallel_bestmodels_all.rds")
message("All done. Combined results saved to predictions_parallel_bestmodels_all.rds")

#--------------------------
# Plot (if final_preds not empty)
#--------------------------
if (nrow(final_preds) > 0) {
  ggplot(final_preds, aes(x = date_predicted)) +
    geom_line(aes(y = actual_y, color = "Actual")) +
    geom_line(aes(y = predicted_y, color = "Predicted")) +
    facet_wrap(~time_horizon, scales = "free_y") +
    labs(title = "Predicted vs Actual (Best Models per Horizon)",
         x = "Date", y = "CPI Value", color = "Legend") +
    theme_minimal()
}
