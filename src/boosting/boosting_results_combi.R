library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)
library(tidyverse)

all_test_preds <- "final_pred_results"


rds_files <- list.files(all_test_preds, pattern = "\\.rds$", full.names = TRUE)

all_preds <- rds_files %>%
  map(readRDS) %>%      # read each .rds file
  bind_rows()           # combine into one big tibble


horizon1 <- all_preds %>%
  filter(time_horizon == 1) %>%
  group_by(dataset) %>%
  mutate(rmse = (mean((predicted_y - actual_y)^2))^0.5) %>%
  mutate(mae = mean(abs(predicted_y - actual_y)))


horizon3 <- all_preds %>%
  filter(time_horizon == 3)%>%
  group_by(dataset) %>%
  slice_head(n=98) %>%
  mutate(rmse = (mean((predicted_y - actual_y)^2))^0.5) %>%
  mutate(mae = mean(abs(predicted_y - actual_y)))


horizon6 <- all_preds %>%
  filter(time_horizon == 6)%>%
  group_by(dataset) %>%
  slice_head(n=95) %>%
  mutate(rmse = (mean((predicted_y - actual_y)^2))^0.5) %>%
  mutate(mae = mean(abs(predicted_y - actual_y)))


horizon12 <- all_preds %>%
  filter(time_horizon == 12)%>%
  group_by(dataset) %>%
  slice_head(n=89) %>%
  mutate(rmse = (mean((predicted_y - actual_y)^2))^0.5) %>%
  mutate(mae = mean(abs(predicted_y - actual_y)))



h1_best <- horizon1 %>%
  group_by(dataset) %>%
  mutate(mean_rmse = mean(rmse, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(mean_rmse == min(mean_rmse)) %>%
  select(dataset, date_predicted, predicted_y, mean_rmse, mae)

h3_best <- horizon3 %>%
  group_by(dataset) %>%
  mutate(mean_rmse = mean(rmse, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(mean_rmse == min(mean_rmse)) %>%
  select(dataset, date_predicted, predicted_y, mean_rmse, mae)

h6_best <- horizon6 %>%
  group_by(dataset) %>%
  mutate(mean_rmse = mean(rmse, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(mean_rmse == min(mean_rmse)) %>%
  select(dataset, date_predicted, predicted_y, mean_rmse, mae)


h12_best <- horizon12 %>%
  group_by(dataset) %>%
  mutate(mean_rmse = mean(rmse, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(mean_rmse == min(mean_rmse)) %>%
  select(dataset, date_predicted, predicted_y, mean_rmse, mae)

saveRDS(list(
  h1_boost = horizon1,
  h3_boost = horizon3,
  h6_boost = horizon6,
  h12_boost = horizon12
), "boost_all.rds")

saveRDS(list(
  h1_best_boost = h1_best,
  h3_best_boost = h3_best,
  h6_best_boost = h6_best,
  h12_best_boost = h12_best
), "boost_best.rds")


library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)  # for combining plots

# Function to generate the plot for any horizon
make_horizon_plot <- function(horizon_data, horizon_label) {
  
  ggplot_dataset_level <- horizon_data %>%
    group_by(dataset) %>%
    select(dataset, date_predicted, predicted_y) %>%
    mutate(y = predicted_y) %>%
    filter(!str_detect(dataset, "_F_"))
  
  ggplot_dataset_y_actual <- horizon_data %>%
    select(date_predicted, actual_y) %>%
    mutate(dataset = "y_real") %>%
    mutate(y = actual_y) %>%
    unique()
  
  ggplot_data <- rbind(ggplot_dataset_level, ggplot_dataset_y_actual)
  
  ggplot(ggplot_data, aes(x = date_predicted, y = y, color = dataset, group = dataset)) +
    geom_line() +
    scale_color_manual(
      values = c("y_real" = "black", "Z_Ht" = "blue", "Z_X_MARX" = "light blue", "Z_X" = "pink")
    ) +
    labs(
      title = paste("Forecasts Horizon", horizon_label),
      x = "Date",
      y = "Predicted Y",
      color = "Dataset"
    ) +
    theme_minimal()
}

# Assuming you have horizon1, horizon3, horizon6, horizon12 data frames:
plot_h1  <- make_horizon_plot(horizon1,  "1")
plot_h3  <- make_horizon_plot(horizon3,  "3")
plot_h6  <- make_horizon_plot(horizon6,  "6")
plot_h12 <- make_horizon_plot(horizon12, "12")

# Combine the 4 plots into one figure (2x2 grid)
combined_plot <- (plot_h1 | plot_h3) /
  (plot_h6 | plot_h12)

# Display combined plot
combined_plot

  