library(gbm)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tibble)
library(parallel)
library(stringr)

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

# We do not touch this for forecasting
y_test <- y_var %>%
  slice((n()-100)+1: n())

# Helper Functions
data_adjusted_time_horizon <- function(data_matrix, time_horizon){
  data_matrix %>%
    mutate(target_var = dplyr::lead(target_var, time_horizon),
           sasdate = dplyr::lead(sasdate, time_horizon)) %>%
    na.omit()
}


#Finding a ballpark for n.trees for all datasets and horizon using OOB 
# parameters_ballpark <- list()
# for (dataset_ in seq_along(all_Z_matrices)){
#   dataset <- all_Z_matrices[[dataset_]]  
#   name_dataset <- names(all_Z_matrices)[dataset_]
#   for (horizon in c(1,3,6,12)){
#     data_matrix_boosting <- dataset %>%
#       mutate(sasdate = ymd(sasdate)) %>%
#       left_join(y_var, by = "sasdate") %>%
#       data_adjusted_time_horizon(horizon)
# 
#     train <- data_matrix_boosting %>%
#       na.omit() %>%
#       slice(1:(n() - 100)) %>%
#       select(-sasdate)
#     boostfit = gbm(target_var~.,data=train,distribution='gaussian',bag.fraction = .5,
#                    interaction.depth=2,n.trees=10000,shrinkage=.1)
#     best = gbm.perf(boostfit, method="OOB")
#     parameters_ballpark[[paste0(name_dataset, "_h", horizon)]] <- list(
#       model = boostfit,
#       best_trees = best
#   }
# }


#------------ FOR datasets Z_F_MARX, Z_F_MARX_Level ---------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

rolling_window_calculate_PCA_train_MARX <- function(overall_train, var_threshold = 0.7) {
  #z_x_cols <- grep("_Z_X$", names(overall_train), value = TRUE)
  z_x_data <- overall_train %>%
    select(-matches("MARX")) %>%
    select(-sasdate)
  
  pca_res <- prcomp(z_x_data, center = TRUE, scale. = TRUE)
  
  explained_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_pc <- which(explained_var >= var_threshold)[1]
  
  pc_scores <- as.data.frame(pca_res$x[, 1:n_pc, drop = FALSE])
  names(pc_scores) <- paste0("PC", 1:n_pc)
  
  overall_train_pca <- overall_train %>%   
    bind_cols(pc_scores)            # add PCA scores
  
  return(list(
    train_pca = overall_train_pca,
    pca_loadings = pca_res$rotation[, 1:n_pc, drop = FALSE],
    n_pc = n_pc,
    explained_variance = explained_var[n_pc]
  ))
}


predict_with_set_parameter_rolling_window_pca <- function(data_matrix_boosting, parameters_boost, time_horizon, dataset_name, test_start_date = as.Date("2012-03-01"), y) {
  
  print(data_matrix_boosting)
  data_matrix_boosting <- data_matrix_boosting %>%
    mutate(sasdate = ymd(sasdate)) %>%
    right_join(y, by = "sasdate") %>%
    data_adjusted_time_horizon(time_horizon) 
  
  print("HAHAHA")
  
  
  test_set <- data_matrix_boosting %>% filter(sasdate >= test_start_date) %>% as_tibble()
  train_set <- data_matrix_boosting %>% filter(sasdate < test_start_date) %>% as_tibble()
  
  if (nrow(test_set) == 0) stop("No test rows found for given test_start_date")
  
  preds_df <- tibble(
    dataset = character(),
    time_horizon = numeric(),
    date_predicted = as.Date(character()),
    predicted_y = numeric(),
    actual_y = numeric(), 
    parameter = numeric()
  )
  
  train_set_full <- train_set
  
  
  for (k in seq_len(nrow(test_set))) {
    
    date_y_t <- test_set$sasdate[k]
    print(date_y_t)
    
    #print(train_set_full)
    pca_res <- rolling_window_calculate_PCA_train_MARX(train_set_full)
    
    #print("FDJFLJSFF")
    
    train_pca <- pca_res$train_pca
    
    boostfit <- gbm(
      target_var ~ .,
      data = train_pca %>% select(-sasdate),  # remove sasdate if exists
      distribution = "gaussian",
      bag.fraction = 0.5,
      interaction.depth = 2,
      n.trees = parameters_boost,
      shrinkage = 0.1
    )
    
    test_row <- test_set[k, , drop = FALSE]
    
    #z_x_cols <- test_row %>%
    #select(-matches("MARX"))
    
    test_z_x <- test_row %>%
      select(-matches("MARX")) %>%
      select(-sasdate)
    
    test_pc <- as.matrix(scale(test_z_x, center = TRUE, scale = TRUE)) %*% pca_res$pca_loadings
    
    print(ncol(test_pc))
    
    test_row_pca <- test_row %>% bind_cols(as_tibble(test_pc))
    names(test_row_pca)[(ncol(test_row_pca) - ncol(test_pc) + 1):ncol(test_row_pca)] <- paste0("PC", seq_len(ncol(test_pc)))
    
    pred_y <- predict(boostfit, newdata = test_row_pca %>% select(-sasdate), n.trees = parameters_boost)
    print(ncol(test_row_pca))
    actual_y <- test_row$target_var
    
    print(pred_y)
    print(actual_y)
    
    preds_df <- add_row(
      preds_df,
      dataset = dataset_name,
      time_horizon = time_horizon,
      date_predicted = date_y_t,
      predicted_y = pred_y,
      actual_y = actual_y, 
      parameter = parameters_boost
    )
    
    save_path <- paste0("predictions_set_para", dataset_name, "_h", time_horizon, parameters_boost, ".rds")
    saveRDS(preds_df, file = save_path)
    
    if (nrow(train_set_full) > 1) {
      train_set_full <- bind_rows(tail(train_set_full, nrow(train_set_full) - 1), test_row)
    } else {
      train_set_full <- test_row
    }
    
  }
  
  return(preds_df)
}

#Parallelized rolling prediction
combinations <- expand.grid(
  datasets = c("Z_F_X_MARX_naked", "Z_F_X_MARX_Level_naked"), 
  horizon = c(1,3,6,12),
  params = c(300, 500, 900)
)

# Function wrapper to run a single combination
run_rolling_combination <- function(combo_row) {
  horizon <- combo_row$horizon
  params  <- combo_row$params
  datasets <- combo_row$datasets
  
  print(paste0("Running horizon=", horizon, ", params=", params))
  
  preds <- predict_with_set_parameter_rolling_window_pca(
    data_matrix_boosting = all_Z_matrices[[datasets]],
    parameters_boost = params,
    time_horizon = horizon,
    dataset_name = datasets,
    test_start_date = as.Date("2012-03-01"), 
    y = y_train
  )
  
  return(preds)
}

# Parallel execution
n_cores <- detectCores() - 1  # leave 1 core free
results <- mclapply(
  seq_len(nrow(combinations)),
  function(idx) run_rolling_combination(combinations[idx, ]),
  mc.cores = n_cores
)

all_preds <- bind_rows(results)

#Find best parameter (lowest MSE) for each prediction horizon, each dataset
best_para_on_preds <- all_preds %>%
  group_by(dataset, time_horizon, parameter) %>%
  filter(row_number() <= 60 - first(time_horizon) + 1) %>%
  mutate(
    rmse = sqrt(mean((predicted_y - actual_y)^2, na.rm = TRUE)),
    mae  = mean(abs(predicted_y - actual_y), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(dataset, time_horizon) %>%
  slice_min(rmse) %>%
  select(dataset, time_horizon, parameter, rmse, mae) %>%
  ungroup() %>%
  unique()

# Running Prediction for test set for each horizon and dataset using parameter chosen above 
for (j in 1:nrow(best_para_on_preds)){
  row <- best_para_on_preds[j,]
  dataset <- row$dataset
  print(dataset)
  parameter <- row$parameter
  print(parameter)
  
  horizon <- row$time_horizon
  print(horizon)
  
  print(all_Z_matrices[[dataset]])

  preds <- predict_with_set_parameter_rolling_window_pca(
    data_matrix_boosting = all_Z_matrices[[dataset]],
    parameters_boost = parameter,
    time_horizon = horizon,
    dataset_name = paste0(dataset, "test_set"),
    test_start_date = as.Date("2017-03-01"), 
    y = y_var
  )
}

#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------



#---------------------- FOR datasets Z_X, Z_Ht, Z_X_MARX ----------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
predict_with_set_parameter <- function(data_matrix_boosting, parameters_boost, time_horizon, dataset_name, test_start_date = as.Date("2012-03-01"), y) {

    data_matrix_boosting <- data_matrix_boosting %>%
    mutate(sasdate = ymd(sasdate)) %>%
    right_join(y, by = "sasdate") %>%
    data_adjusted_time_horizon(time_horizon)
  
  test_set <- data_matrix_boosting %>% filter(sasdate >= test_start_date) %>% as_tibble()
  train_set <- data_matrix_boosting %>% filter(sasdate < test_start_date) %>% as_tibble()
  
  preds_df <- tibble(
    dataset = character(),
    time_horizon = numeric(),
    date_predicted = as.Date(character()),
    predicted_y = numeric(),
    actual_y = numeric(),
    parameter = numeric(),
  )
  
  initial_train_size <- nrow(train_set)
  train_set <- train_set %>%
    select(-sasdate)
  test_set_ <- test_set %>%
    select(-sasdate)
  
  # iterate through test rows and do rolling-origin updates
  for (k in seq_len(nrow(test_set))) {
    date_y_t <- test_set$sasdate[k]
    
    # get best model on current train_set
    boostfit = gbm(target_var~.,data=train_set,distribution='gaussian',bag.fraction = .5,
                   interaction.depth=2,n.trees=parameters_boost,shrinkage=.01)
    
    
    
    # predict for current test row
    newdata_row <- test_set[k, setdiff(names(test_set), "target_var"), drop = FALSE] %>% as_tibble()
    pred_y <- predict(boostfit, newdata = newdata_row, n.trees = parameters_boost)
    actual_y <- test_set$target_var[k]
    
    print(pred_y)
    print(actual_y)
    print(date_y_t)
    
    # add to preds_df
    preds_df <- add_row(preds_df,
                        dataset = dataset_name,
                        time_horizon = time_horizon,
                        date_predicted = date_y_t,
                        predicted_y = pred_y,
                        actual_y = actual_y,
                        parameter = parameters_boost
    )
    save_path <- paste0("predictions_set_para", dataset_name, "_h", time_horizon, parameters_boost, ".rds")
    saveRDS(preds_df, file = save_path)
    
    # Expand training set for rolling-origin safely
    
    if (nrow(train_set) > 1) {
      train_set <- bind_rows(tail(train_set, nrow(train_set) - 1) %>% as_tibble(), test_set_[k, , drop = FALSE] %>% as_tibble())
    } else {
      train_set <- test_set[k, , drop = FALSE] %>% as_tibble()
    }
  }
  
  return(preds_df)
}

#Parallelized rolling prediction
combinations <- expand.grid(
  datasets = c("Z_Ht", "Z_X_MARX", "Z_X"), 
  horizon = c(1,3,6,12),
  params = c(3, 5, 9)
)

# Function wrapper to run a single combination
run_rolling_combination <- function(combo_row) {
  horizon <- combo_row$horizon
  params  <- combo_row$params
  datasets <- combo_row$datasets
  
  print(paste0("Running horizon=", horizon, ", params=", params))
  
  preds <- predict_with_set_parameter(
    data_matrix_boosting = all_Z_matrices[[datasets]],
    parameters_boost = params,
    time_horizon = horizon,
    dataset_name = datasets,
    test_start_date = as.Date("2012-03-01"), 
    y = y_train
  )
  
  return(preds)
}

# Parallel execution
n_cores <- detectCores() - 1  # leave 1 core free
results <- mclapply(
  seq_len(nrow(combinations)),
  function(idx) run_rolling_combination(combinations[idx, ]),
  mc.cores = n_cores
)
all_preds <- bind_rows(results)
#saveRDS(all_preds, "Z_Ht_final_overall.rds")

#Find best parameter (lowest MSE) for each prediction horizon, each dataset
best_para_on_preds <- all_preds %>%
  group_by(dataset, time_horizon, parameter) %>%
  filter(row_number() <= 60 - first(time_horizon) + 1) %>%
  mutate(
    rmse = sqrt(mean((predicted_y - actual_y)^2, na.rm = TRUE)),
    mae  = mean(abs(predicted_y - actual_y), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(dataset, time_horizon) %>%
  slice_min(rmse) %>%
  select(dataset, time_horizon, parameter, rmse, mae) %>%
  ungroup() %>%
  unique()

# Running Prediction for test set for each horizon and dataset using parameter chosen above 
for (j in 1:nrow(best_para_on_preds)){
  row <- best_para_on_preds[j,]
  dataset <- best_para_on_preds$dataset
  parameter <- best_para_on_preds$parameter
  horizon <- best_para_on_preds$time_horizon
  
  preds <- predict_with_set_parameter(
    data_matrix_boosting = all_Z_matrices[[dataset]],
    parameters_boost = params,
    time_horizon = horizon,
    dataset_name = paste0(dataset, "test_set"),
    test_start_date = as.Date("2017-03-01"), 
    y = y_var
  )
}

# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------

# ---------------------------- FOR datasets F_X, F_Level ------------------------------------------------
# -------------------------------------------------------------------------------------------------------
rolling_window_calculate_PCA_train_F_X <- function(overall_train, var_threshold = 0.7) {
  
  z_x_cols <- grep("_Z_X$", names(overall_train), value = TRUE)
  z_x_data <- overall_train %>% select(all_of(z_x_cols))
  
  pca_res <- prcomp(z_x_data, center = TRUE, scale. = TRUE)
  
  explained_var <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
  n_pc <- which(explained_var >= var_threshold)[1]
  
  pc_scores <- as.data.frame(pca_res$x[, 1:n_pc, drop = FALSE])
  names(pc_scores) <- paste0("PC", 1:n_pc)
  
  overall_train_pca <- overall_train %>%
    select(-all_of(z_x_cols)) %>%   # drop original Z_X columns
    bind_cols(pc_scores)            # add PCA scores
  
  return(list(
    train_pca = overall_train_pca,
    pca_loadings = pca_res$rotation[, 1:n_pc, drop = FALSE],
    n_pc = n_pc,
    explained_variance = explained_var[n_pc]
  ))
}



predict_with_set_parameter_rolling_window_pca_ <- function(data_matrix_boosting, parameters_boost, time_horizon, dataset_name, test_start_date = as.Date("2012-03-01"), y) {
  data_matrix_boosting <- data_matrix_boosting %>%
    mutate(sasdate = ymd(sasdate)) %>%
    right_join(y, by = "sasdate") %>%
    data_adjusted_time_horizon(time_horizon)
  
  full_matrix <- all_Z_matrices[["Z_X"]] %>%
    rename_with(~paste0(.x, "_Z_X"), -sasdate)
  
  overall_train <- full_matrix %>%
    right_join(data_matrix_boosting, by = "sasdate", suffix = c(".Z_X", ".train")) %>%
    na.omit()
  
  #print(overall_train)
  
  
  test_set <- overall_train %>% filter(sasdate >= test_start_date) %>% as_tibble()
  train_set <- overall_train %>% filter(sasdate < test_start_date) %>% as_tibble()
  
  if (nrow(test_set) == 0) stop("No test rows found for given test_start_date")
  
  preds_df <- tibble(
    dataset = character(),
    time_horizon = numeric(),
    date_predicted = as.Date(character()),
    predicted_y = numeric(),
    actual_y = numeric()
  )
  
  train_set_full <- train_set
  
  #print("HAHAHAHHA")
  
  for (k in seq_len(nrow(test_set))) {
    
    date_y_t <- test_set$sasdate[k]
    print(date_y_t)
    
    pca_res <- rolling_window_calculate_PCA_train_F_X(train_set_full)
    
    train_pca <- pca_res$train_pca
    
    boostfit <- gbm(
      target_var ~ .,
      data = train_pca %>% select(-sasdate),  # remove sasdate if exists
      distribution = "gaussian",
      bag.fraction = 0.5,
      interaction.depth = 2,
      n.trees = parameters_boost,
      shrinkage = 0.01
    )
    
    test_row <- test_set[k, , drop = FALSE]
    
    z_x_cols <- grep("_Z_X$", names(test_set), value = TRUE)
    test_z_x <- test_row %>% select(all_of(z_x_cols))
    
    test_pc <- as.matrix(scale(test_z_x, center = TRUE, scale = TRUE)) %*% pca_res$pca_loadings
    
    print(ncol(test_pc))
    
    test_row_pca <- test_row %>% select(-all_of(z_x_cols)) %>% bind_cols(as_tibble(test_pc))
    names(test_row_pca)[(ncol(test_row_pca) - ncol(test_pc) + 1):ncol(test_row_pca)] <- paste0("PC", seq_len(ncol(test_pc)))
    
    pred_y <- predict(boostfit, newdata = test_row_pca %>% select(-sasdate), n.trees = parameters_boost)
    actual_y <- test_row$target_var
    
    print(pred_y)
    print(actual_y)
    
    preds_df <- add_row(
      preds_df,
      dataset = dataset_name,
      time_horizon = time_horizon,
      date_predicted = date_y_t,
      predicted_y = pred_y,
      actual_y = actual_y
    )
    
    save_path <- paste0("predictions_set_para", dataset_name, "_h", time_horizon, parameters_boost, ".rds")
    saveRDS(preds_df, file = save_path)
    
    if (nrow(train_set_full) > 1) {
      train_set_full <- bind_rows(tail(train_set_full, nrow(train_set_full) - 1), test_row)
    } else {
      train_set_full <- test_row
    }
    
  }
  
  return(preds_df)
}

#Parallelized rolling prediction
combinations <- expand.grid(
  datasets = c("Z_F_naked", "Z_F_Level_naked"), 
  horizon = c(1,3,6,12),
  params = c(300, 500, 900)
)

# Function wrapper to run a single combination
run_rolling_combination <- function(combo_row) {
  horizon <- combo_row$horizon
  params  <- combo_row$params
  datasets <- combo_row$datasets
  
  print(paste0("Running horizon=", horizon, ", params=", params))
  
  preds <- predict_with_set_parameter_rolling_window_pca_(
    data_matrix_boosting = all_Z_matrices[[datasets]],
    parameters_boost = params,
    time_horizon = horizon,
    dataset_name = datasets,
    test_start_date = as.Date("2012-03-01"), 
    y = y_train
  )
  
  return(preds)
}

# Parallel execution
n_cores <- detectCores() - 1  # leave 1 core free
results <- mclapply(
  seq_len(nrow(combinations)),
  function(idx) run_rolling_combination(combinations[idx, ]),
  mc.cores = n_cores
)
all_preds <- bind_rows(results)


#Find best parameter (lowest MSE) for each prediction horizon, each dataset
best_para_on_preds <- all_preds %>%
  group_by(dataset, time_horizon, parameter) %>%
  filter(row_number() <= 60 - first(time_horizon) + 1) %>%
  mutate(
    rmse = sqrt(mean((predicted_y - actual_y)^2, na.rm = TRUE)),
    mae  = mean(abs(predicted_y - actual_y), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(dataset, time_horizon) %>%
  slice_min(rmse) %>%
  select(dataset, time_horizon, parameter, rmse, mae) %>%
  ungroup() %>%
  unique()

# Running Prediction for test set for each horizon and dataset using parameter chosen above 
for (j in 1:nrow(best_para_on_preds)){
  row <- best_para_on_preds[j,]
  dataset <- best_para_on_preds$dataset
  parameter <- best_para_on_preds$parameter
  horizon <- best_para_on_preds$time_horizon
  
  preds <- predict_with_set_parameter_rolling_window_pca_(
    data_matrix_boosting = all_Z_matrices[[dataset]],
    parameters_boost = params,
    time_horizon = horizon,
    dataset_name = paste0(dataset, "test_set"),
    test_start_date = as.Date("2017-03-01"), 
    y = y_var
  )
}
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------


