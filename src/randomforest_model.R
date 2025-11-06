# load libraries
library(readr)
library(dplyr)
library(tidyverse)
library(randomForest)

# set seed for reproducibility
set.seed(4308)

# --------- BENCHMARK USING TRANSFORMED PREDICTORS ---------

# load the X values
predictors_transformed <- read_csv("../data/predictors_transformed.csv")
X <- predictors_transformed

# load the y values
cpi_target_full <- read_csv("../data/cpi_target_full.csv")
y <- cpi_target_full[, -2]  # assuming column 2 is dropped

runrf <- function(X, y, lag = 1) {
  # combine both X and y using left join on sasdate
  full_data <- X %>% left_join(y, by = "sasdate") %>% drop_na()
  
  # drop date column and convert to numeric matrix
  data_mat <- as.matrix(full_data[, -1])
  
  # create 4 lags + forecast horizon shift (=lag)
  temp <- embed(data_mat, 4 + lag)
  
  # define number of predictors (including target)
  nvar <- ncol(data_mat)
  
  # y variable (inflation rate target = last column before embedding)
  y_aligned <- temp[, nvar]
  
  # lagged predictors (drop future columns)
  X_lagged <- temp[, -c(1:(nvar * lag))]
  
  # fit the random forest on default settings
  model <- randomForest(X_lagged, y_aligned, importance = TRUE)
  
  # generate forecast (using the last available observation)
  if (nrow(X_lagged) > 0) {
    pred <- tail(predict(model, X_lagged), 1)
  } else {
    pred <- NA
  }
  
  # return the estimated model and h-step forecast
  return(list("model" = model, "pred" = pred))
}

rf_rolling_window <- function(X, y, nprev, lag = 1) {
  # combine into one full dataset
  full_data <- X %>% 
    left_join(y, by = "sasdate") %>% 
    drop_na()
  
  # storage for outputs
  save_importance <- vector("list", nprev)
  save_pred <- matrix(NA, nprev, 1)
  
  # rolling window loop
  for (i in nprev:1) {
    # define the rolling estimation window
    Y_window <- full_data[(1 + nprev - i):(nrow(full_data) - i), ]
    
    # call the RF forecasting function
    rf_result <- runrf(
      X = Y_window[, -ncol(Y_window)],  # all columns except target
      y = Y_window[, c("sasdate", names(Y_window)[ncol(Y_window)])],  # sasdate + target
      lag = lag
    )
    
    # store predictions and variable importance
    save_pred[(1 + nprev - i), ] <- rf_result$pred
    save_importance[[i]] <- importance(rf_result$model)
    
    cat("Iteration", (1 + nprev - i), "completed. \n")
  }
  
  # evaluation data
  real <- full_data[[ncol(full_data)]]
  time_axis <- full_data$sasdate
  pred_padded <- c(rep(NA, length(real) - nprev), save_pred)
  
  # compute forecast errors
  rmse <- sqrt(mean((tail(real, nprev) - save_pred)^2))
  mae <- mean(abs(tail(real, nprev) - save_pred))
  errors <- c("rmse" = rmse, "mae" = mae)
  
  # return everything cleanly
  return(list(
    pred = data.frame(sasdate = time_axis,
                      actual = real,
                      forecast = pred_padded),
    errors = errors,
    save_importance = save_importance
  ))
}

# plotting function
plot_rf_forecast <- function(rf_results, lag_val) {
  # extract variables
  time_axis <- as.Date(rf_results$pred$sasdate)
  real <- rf_results$pred$actual
  pred_padded <- rf_results$pred$forecast
  rmse_val <- round(rf_results$errors["rmse"], 4)
  mae_val  <- round(rf_results$errors["mae"], 4)
  
  # titles
  main_title <- sprintf("RF Rolling Forecasts (Lag = %d Months)", lag_val)
  sub_title  <- sprintf("RMSE = %.4f | MAE = %.4f", rmse_val, mae_val)
  
  # base R plot
  plot(time_axis, real, type = "l",
       main = main_title, sub = sub_title,
       ylab = "Inflation", xlab = "Time",
       col = "black", lwd = 1.2, xaxt = "n")
  
  # add forecast line
  lines(time_axis, pred_padded, col = "red", lwd = 1.2)
  
  # add custom date axis
  axis.Date(side = 1, at = seq(min(time_axis), max(time_axis), by = "1 year"),
            format = "%Y")
  
  # legend and grid
  legend("topleft", legend = c("Actual", "Forecast"),
         col = c("black", "red"), lty = 1, bty = "n", horiz = TRUE)
  grid(col = "gray80")
}


# run the random forest with rolling window
rf_results <- rf_rolling_window(
  X = predictors_transformed,
  y = cpi_target_full %>% select(sasdate, CPI_t),
  nprev = 66,  # 66 out-of-sample forecasts (starting from 2020 Jan)
  lag = 1      # 1-month ahead
)

# plot graph
plot_rf_forecast(rf_results, 1)

# view results
print(rf_results$errors)


# --------- SAVE RESULTS ---------
saveRDS(rf_results, file = "../data/rf_results/predictors_transformed_lag1.rds")


# --------- LOOP THROUGH FOR THE DIFFERENT LAGS ---------
# run the random forest with rolling window across the different forecast horizons
for (lag_val in c(1, 3, 6, 12)) {
  cat("Starting Random Forest Rolling Forecast for Lag =", lag_val, "\n")
  cat("🕒 Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  # run rolling window random forest
  rf_results <- rf_rolling_window(
    X = predictors_transformed,
    y = cpi_target_full %>% select(sasdate, CPI_t),
    nprev = 66,    # 66 out-of-sample forecasts (starting from 2020 Jan)
    lag = lag_val  # forecast horizon
    )
  
  cat("Model Completed for Lag =", lag_val, "\n")
  
  # save the full rf_results object
  result_filename <- sprintf("../data/rf_results/predictors_transformed_lag%d.rds", lag_val)
  saveRDS(rf_results, file = result_filename)
  
  # save plot
  png_filename <- sprintf("../figures/rf_results/predictors_transformed_lag%d_plot.png", lag_val)
  png(png_filename, width = 800, height = 500)
  plot_rf_forecast(rf_results, lag_val)
  dev.off()
  
  # print key metrics
  if (!is.null(rf_results$errors)) {
    print(rf_results$errors)
  }

  cat("Saved Results to:", result_filename, "\n")
  
  cat("Finished Lag =", lag_val, "\n")
  cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
}
  
  




# unpack all the Z matrices from feature engineering
Z_objects <- readRDS("all_Z_matrices.rds")
list2env(Z_objects, envir = .GlobalEnv)
Z_matrices <- names(Z_objects)