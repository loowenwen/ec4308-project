library(tidyverse)
library(lubridate)
library(timetk)
library(randomForest)
library(caret)
library(scales)
set.seed(4308)

# ----------------------------------------------------------
# 1. Load and inspect raw FRED-MD data
# ----------------------------------------------------------
# The FRED-MD monthly dataset usually has the first column as dates
# and remaining columns as macro variables.
# CPIAUCSL = Consumer Price Index (All Urban Consumers)
# ----------------------------------------------------------
fred_raw <- read_csv("../data/2025-09-MD.csv")

glimpse(fred_raw[1:10])  # quick peek

# convert date column (usually named 'sasdate' or 'date')
if ("sasdate" %in% names(fred_raw)) {
  fred_raw <- fred_raw %>% rename(date = sasdate)
}
# remove first row
fred_raw <- fred_raw[-1,]
fred_raw$date <- as.Date(fred_raw$date, format = "%m/%d/%Y")

# ----------------------------------------------------------
# 2. Construct target variable: Monthly CPI inflation (% change)
# ----------------------------------------------------------
# CPIAUCSL is usually in levels, so we take log-difference * 1200
# (to get annualized monthly %)
# ----------------------------------------------------------
fred_raw <- fred_raw %>%
  mutate(inflation = 1200 * (log(CPIAUCSL) - lag(log(CPIAUCSL)))) %>%
  drop_na(inflation)

# ----------------------------------------------------------
# 3. Clean predictors
# ----------------------------------------------------------
# Remove CPIAUCSL (used to create inflation)
# Remove any non-numeric or constant columns
# ----------------------------------------------------------
fred_clean <- fred_raw %>%
  select(date, inflation, where(is.numeric)) %>%
  select(-CPIAUCSL) %>%
  mutate(across(-c(date, inflation), scale)) %>%
  drop_na()

# ----------------------------------------------------------
# 4. Feature engineering: add inflation lags (1:12 months)
# ----------------------------------------------------------
fred_lagged <- fred_clean %>%
  tk_augment_lags(.value = inflation, .lags = 1:12) %>%
  drop_na()

# ----------------------------------------------------------
# 5. Train-test split (train: pre-2020, test: 2020+)
# ----------------------------------------------------------
train <- fred_lagged %>% filter(date < "2020-01-01")
test  <- fred_lagged %>% filter(date >= "2020-01-01")

X_train <- train %>% select(-c(date, inflation))
y_train <- train$inflation
X_test  <- test %>% select(-c(date, inflation))
y_test  <- test$inflation

cat("Training sample:", nrow(train), "obs\n")
cat("Testing sample:", nrow(test), "obs\n")

# ----------------------------------------------------------
# 6. Train Random Forest
# ----------------------------------------------------------
rf_model <- randomForest(
  x = X_train,
  y = y_train,
  ntree = 500,
  mtry = floor(sqrt(ncol(X_train))),
  importance = TRUE
)

print(rf_model)

# ----------------------------------------------------------
# 7. Forecast & evaluate
# ----------------------------------------------------------
pred_rf <- predict(rf_model, X_test)

results <- tibble(
  date = test$date,
  actual = y_test,
  forecast = pred_rf,
  error = forecast - actual
)

rmse_rf <- sqrt(mean(results$error^2))
mae_rf  <- mean(abs(results$error))
direction_acc <- mean(sign(results$forecast) == sign(results$actual))

cat("\nRandom Forest performance:")
cat("\nRMSE:", round(rmse_rf, 4))
cat("\nMAE :", round(mae_rf, 4))
cat("\nDirectional Accuracy:", round(direction_acc, 3), "\n")

# ----------------------------------------------------------
# 8. Visualization
# ----------------------------------------------------------
results %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = actual, colour = "Actual Inflation"), linewidth = 0.9) +
  geom_line(aes(y = forecast, colour = "RF Forecast"), linewidth = 0.9, linetype = "dashed") +
  geom_vline(xintercept = as.Date("2020-01-01"), linetype = "dotted", color = "black") +
  annotate("text", x = as.Date("2020-06-01"), y = max(results$actual, na.rm=TRUE)*0.9,
           label = "COVID-19 period", size = 3.5, color = "gray30") +
  scale_colour_manual(values = c("Actual Inflation"="black","RF Forecast"="firebrick")) +
  labs(title = "Random Forest Forecast vs Actual Inflation",
       subtitle = "Training: 1980–2019 | Testing: 2020–2024",
       y = "Inflation (annualized, %)", x = "", colour = "") +
  theme_minimal(base_size = 13)

# ----------------------------------------------------------
# 9. Variable importance
# ----------------------------------------------------------
importance_df <- as.data.frame(importance(rf_model)) %>%
  rownames_to_column("Variable") %>%
  arrange(desc(IncNodePurity)) %>%
  slice(1:15)

ggplot(importance_df, aes(x = reorder(Variable, IncNodePurity),
                          y = IncNodePurity)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 15 Important Predictors (Random Forest)",
       x = "", y = "Increase in Node Purity") +
  theme_minimal()