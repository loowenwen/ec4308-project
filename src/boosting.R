library(gbm)
library(lubridate)
library(ggplot2)

fred_data_clean <- read.csv("../data/fred_data_clean.csv")

y_var = fred_data_clean %>%
  select(CPIAUCSL, sasdate) %>%
  mutate(sasdate = ymd(sasdate)) %>%
  mutate(target_var = CPIAUCSL) %>%
  select(-CPIAUCSL)

# Datasets to run Gradient Boosting on 
datasets <- c("Z_12")

# Parameter to choose for depth of Tree
parameters_boost = c(15)

# Helper Function to adjust values for Time Horizon 
data_adjusted_time_horizon <- function(data_matrix, time_horizon){
  data_matrix_new <- data_matrix%>%
    mutate(target_var = dplyr::lead(target_var, time_horizon-1)) %>%
    mutate(sasdate = dplyr::lead(sasdate, time_horizon-1)) %>%
    na.omit()
  return(data_matrix_new)
}

# Helper function that loops through validation set to find best parameter (depth of Tree)
boosting_for_one_parameter <- function(train_boosting, val_boosting, parameter){
  mse_matrix = matrix(0, nrow(val_boosting), 1)
  for (i in 1:nrow(val_boosting)){
    boostfit = gbm(target_var ~., data = train_boosting, distribution='gaussian', bag.fraction=.5,
                   interaction.depth=2, n.trees=parameter, shrinkage=.1)
    #print(i)
    value_predicted = predict(boostfit, newdata = val_boosting[i,], n.trees=parameter)
    print(value_predicted)
    mse = (value_predicted - val_boosting$target_var[i])^2
    #print(mse)
    mse_matrix[i] = mse 
    train_boosting = rbind(train_boosting[2:nrow(train_boosting),], val_boosting[i,])
    #print(nrow(train_boosting))
    
  }
  return(list(mse = sqrt(mean(mse_matrix, na.rm = TRUE)), fitted_model = boostfit))
}

# Table to store all models with parameters and time horizon 
models_tbl_overall <- tibble(
  dataset_name = character(),
  model = list(),
  mse = numeric(),
  parameters = numeric(),
  y_t = Date(), 
  time_horizon = numeric()
)

# Main Gradient Boosting Algorithm 
gradient_boosting <- function(data_matrix_boost, parameters_boost, time_horizon, dataset_name, date){
  boosting_matrix <- data_matrix_boost %>%
    left_join(y_var, by = "sasdate")
  boosting_matrix <- data_adjusted_time_horizon(boosting_matrix, time_horizon)
  boosting_matrix <- boosting_matrix %>%
    select(-sasdate)
  train_boosting = boosting_matrix[1:173,]
  val_boosting = boosting_matrix[174:nrow(boosting_matrix),]
  for (j in parameters_boost){
    res = boosting_for_one_parameter(train_boosting, val_boosting, j)
    models_tbl_overall <<- add_row(models_tbl_overall, 
                           dataset_name = dataset_name, 
                           model = list(res$fitted_model), 
                           mse = res$mse, 
                           parameters = j,
                           y_t = date,
                           time_horizon = time_horizon 
                           )
  }
  
  return(models_tbl_overall)
}

#for (i in datasets){
  #dataset <- get(i)
  #matrix = gradient_boosting(dataset, parameters_boost, 20, i)
#}

predict_test_set_per_time_horizon <- function(data_matrix_boosting, parameters_boost, time_horizon, dataset_name){
  test_set <- data_matrix_boosting %>%
    filter(sasdate >= as.Date("2020-01-01")) 
  train_set <- data_matrix_boosting %>%
    filter(sasdate < as.Date("2020-01-01")) 
  for (k in 1:(nrow(test_set) - time_horizon)){
    print(k)
    gradient_boosting(train_set, parameters_boost, time_horizon, dataset_name, test_set$sasdate[k])
    train_set <- rbind(train_set[(1+time_horizon):nrow(train_set), ], test_set[k+time_horizon, ])
  }
}

predict_test_set_per_time_horizon(Z_12,parameters_boost,1,"Z_12")

# Find the best model for each time horizon 
best_model_each_time_horizon <- models_tbl_overall %>%
  group_by(time_horizon, y_t) %>%
  mutate(min_mse = min(mse, na.rm = TRUE)) %>%
  filter(mse == min_mse) %>%
  arrange(time_horizon)

# Data Matrix with target Var
data_matrix <- Z_12 %>%
  left_join(y_var, by = "sasdate")

best_model_each_time_horizon$predicted_y <- NA
best_model_each_time_horizon$actual_y <- NA

horizon = 0 
## To predict each horizon
for (row in 1:nrow(best_model_each_time_horizon)){
  row_i <- best_model_each_time_horizon[row,]
  time_horizon = row_i$time_horizon
  if (time_horizon != horizon){
    temp_matrix <- data_matrix %>% right_join(best_model_each_time_horizon, by = c("sasdate" = "y_t"))
    horizon = time_horizon
  }
  
  time_horizon <- row_i$time_horizon 
  model <- row_i$model[[1]]
  parameter <- row_i$parameters
  best_model_each_time_horizon$predicted_y[row] = predict(model, newdata = temp_matrix[row, ], n.trees=parameter)
  best_model_each_time_horizon$actual_y[row] = temp_matrix$target_var[row]
  
}

## Plot

ggplot(best_model_each_time_horizon, aes(x = y_t)) +
  geom_line(aes(y = actual_y, color = "Actual")) +
  geom_line(aes(y = predicted_y, color = "Predicted")) +
  labs(
    title = "Predicted vs Actual over Time",
    x = "Date",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal()