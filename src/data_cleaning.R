###############################################################
# Data Cleaning and Transformation Script
# - Handles missing data
# - Applies FRED-MD transformations (1–7)
# - Forward-fills short gaps (e.g., Apr 2020)
# - Checks for stationarity
# - Splits dataset into training (pre-2020) and testing (post-2020)
###############################################################

# ---------------------- 1. Load Libaries ---------------------
library(zoo)
library(readr)
library(dplyr)
library(lubridate)
library(tseries)
library(tidyverse)
library(urca)
library(ISLR)
library(pls)
library(glmnet)
attach(Credit)

# ---------------------- 2. Import Raw Data -------------------
# Read FRED-MD CSV (update path accordingly)
fred_raw <- read_csv("../data/2025-09-MD.csv")
head(fred_raw)

# Standardize column names and parse dates
fred_data <- fred_raw %>%
  rename(sasdate = 1) %>%
  mutate(sasdate = as.Date(sasdate, format = "%m/%d/%Y"))

# ---------------------- 3. Filter Sample Window ---------------
# Keep data from Jan 1960 onwards
fred_data <- fred_data %>%
  filter(sasdate >= as.Date("1960-01-01"))

# ---------------------- 4. Inspect Missingness ----------------
# Summarize missing values before transformation
na_summary <- data.frame(
  variable = names(fred_data)[-1],
  n_missing = colSums(is.na(fred_data[-1])),
  pct_missing = round(100 * colMeans(is.na(fred_data[-1])), 2)
) %>%
  arrange(desc(pct_missing))

head(na_summary, 15) # preview top 15 variables with most missing data

# ---------------------- 5. Drop Heavily Incomplete Columns ----
# Drop problem series with high NA rates
drop_vars <- c("ACOGNO", "TWEXAFEGSMTHx", "UMCSENTx", "ANDENOx", "VIXCLSx")

fred_data <- fred_data %>%
  select(-all_of(drop_vars))

# ---------------------- 6. Identify Remaining Missing Rows ----
# Identify rows and columns with any remaining NAs
fred_data %>%
  filter(if_any(everything(), is.na))

# ---------------------- 7. Review Specific NA Dates -----------
# Known problematic dates: Apr 2020, Jul 2025, Aug 2025
target_dates <- as.Date(c("2020-04-01", "2025-07-01", "2025-08-01"))

missing_rows <- fred_data %>%
  filter(sasdate %in% target_dates)

missing_map <- missing_rows %>%
  pivot_longer(cols = -sasdate, names_to = "variable", values_to = "value") %>%
  filter(is.na(value)) %>%
  arrange(sasdate)

print(missing_map)

# ---------------------- 8. Handle Incomplete Months -----------
# Drop incomplete future months (Jul & Aug 2025)
fred_clean <- fred_data %>%
  filter(sasdate <= as.Date("2025-06-01"))

# Forward fill short one-month gaps (Apr 2020) BEFORE transformation
fred_clean <- fred_clean %>%
  mutate(
    CP3Mx = na.locf(CP3Mx, na.rm = FALSE),
    COMPAPFFx = na.locf(COMPAPFFx, na.rm = FALSE)
  )

# Quick check that no NAs remain in those target months
fred_clean %>%
  filter(sasdate %in% target_dates) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "num_missing"
  ) %>%
  filter(num_missing > 0)


# ---------------------- 9. Define Transformation Functions ----
# FRED-MD transformation codes (1–7)

# Code 1: no transformation (level)
transform_1 <- function(x) {
  x
}

# Code 2: first difference
transform_2 <- function(x) {
  x - lag(x, n = 1, default = NA)
}

# Code 3: second difference
transform_3 <- function(x) {
  (x - lag(x, 1)) - (lag(x, 1) - lag(x, 2))
}

# Code 4: log transform
transform_4 <- function(x) {
  ifelse(x > 0, log(x), NA) # safer than plain log(x)
}

# Code 5: log first difference
transform_5 <- function(x) {
  log(x) - lag(log(x), 1)
}

# Code 6: log second difference
transform_6 <- function(x) {
  (log(x) - lag(log(x), 1)) - (lag(log(x), 1) - lag(log(x), 2))
}

# Code 7: custom ratio transform (rarely used)
transform_7 <- function(x) {
  (x / (lag(x, 1) - 1)) - (lag(x, 1) / (lag(x, 2) - 1))
}


# ---------------------- 10. Apply Transformations -------------
# Preserve a copy of the cleaned levels prior to differencing
fred_complete_levels <- fred_clean
write_csv(fred_complete_levels, "../data/fred_complete_levels.csv")

# Extract and save transform codes
transform_codes_full <- fred_raw[1, -1] %>%
  unlist(use.names = TRUE) %>%
  as.numeric()
names(transform_codes_full) <- names(fred_raw)[-1]
write_csv(tibble(variable = names(transform_codes_full), transform_code = transform_codes_full), "../data/transform_codes_.csv")

# Align transformation codes with currently retained columns
cols_to_transform <- names(fred_clean)[-1]
transform_codes <- transform_codes_full[cols_to_transform]
if (all(is.na(transform_codes))) warning("All inferred transform codes are NA. Please confirm the raw file contains the transformation row as expected.")

# Force CPI to transform code 5 (log first difference) per instructor
if ("CPIAUCSL" %in% names(transform_codes)) {
  transform_codes[["CPIAUCSL"]] <- 5
  transform_codes_full[["CPIAUCSL"]] <- 5
  write_csv(tibble(variable = names(transform_codes_full), transform_code = transform_codes_full), "../data/transform_codes_.csv")
}

# Initialize transformed data from cleaned levels
fred_data_values <- fred_clean
transforms <- list(transform_1, transform_2, transform_3, transform_4, transform_5, transform_6, transform_7)

# Apply transforms to all columns except CPI (we'll compute CPI target separately)
for (col_name in cols_to_transform) {
  if (col_name == "CPIAUCSL") next
  code <- transform_codes[[col_name]]
  if (is.na(code) || code < 1 || code > length(transforms)) next
  fn <- transforms[[code]]
  fred_data_values[[col_name]] <- fn(as.numeric(fred_data_values[[col_name]]))
}

# Compute CPI target from raw cleaned levels: log first difference (percent)
if ("CPIAUCSL" %in% names(fred_complete_levels)) {
  cpi_raw <- as.numeric(fred_complete_levels$CPIAUCSL)
  cpi_t <- 100 * transform_5(cpi_raw)
} else {
  stop("CPIAUCSL not found in cleaned levels; cannot compute CPI target.")
}

# Determine rows lost to differencing/logging and drop them for all series
loss_map <- c(`1` = 0, `2` = 1, `3` = 2, `4` = 0, `5` = 1, `6` = 2, `7` = 12)
max_loss_preds <- max(loss_map[as.character(transform_codes)], na.rm = TRUE)
max_loss_cpi <- loss_map[as.character(5)]
max_loss <- max(max_loss_preds, max_loss_cpi, na.rm = TRUE)
if (is.finite(max_loss) && max_loss > 0) {
  message("Dropping first ", max_loss, " rows lost to differencing/logging.")
  fred_data_values <- fred_data_values %>% slice(-(1:max_loss))
  cpi_raw <- cpi_raw[-(1:max_loss)]
  cpi_t <- cpi_t[-(1:max_loss)]
}

# Drop any rows with remaining NAs in predictors
fred_data_values <- fred_data_values %>% drop_na()

# Align CPI levels by date (join on sasdate) to ensure exact matching rows
cpi_levels_sub <- fred_complete_levels %>%
  filter(sasdate %in% fred_data_values$sasdate) %>%
  arrange(sasdate)
fred_data_values <- fred_data_values %>% arrange(sasdate)

if (nrow(cpi_levels_sub) != nrow(fred_data_values)) {
  # In case of mismatch, inner-join to be safe
  merged <- inner_join(fred_data_values %>% select(sasdate), fred_complete_levels %>% select(sasdate, CPIAUCSL), by = "sasdate")
  cpi_raw_aligned <- as.numeric(merged$CPIAUCSL)
  cpi_t_aligned <- 100 * transform_5(cpi_raw_aligned)
  cpi_df <- tibble(sasdate = merged$sasdate, CPI_raw = cpi_raw_aligned, CPI_t = cpi_t_aligned)
} else {
  cpi_raw_aligned <- as.numeric(cpi_levels_sub$CPIAUCSL)
  cpi_t_aligned <- 100 * transform_5(cpi_raw_aligned)
  cpi_df <- tibble(sasdate = cpi_levels_sub$sasdate, CPI_raw = cpi_raw_aligned, CPI_t = cpi_t_aligned)
}

# Predictors (exclude CPI) and save
predictors_df <- fred_data_values %>% select(-any_of("CPIAUCSL"))
write_csv(predictors_df, "../data/predictors_transformed.csv")
write_csv(cpi_df, "../data/cpi_target_full.csv")


# ---------------------- 11. Stationarity Checks ---------------
# Augmented Dickey-Fuller (ADF) test without trend
list_of_non_stationary <- list()

# Loop through each variable (excluding date)
for (i in 2:ncol(fred_data_values)) {
  adf_p <- tryCatch(adf.test(fred_data_values[[i]], k = 0)$p.value,
    error = function(e) NA
  )
  if (!is.na(adf_p) && adf_p > 0.05) {
    list_of_non_stationary <- c(list_of_non_stationary, i)
  }
}

# ADF test with trend
list_of_non_stationary_trend <- list()
for (i in 2:ncol(fred_data_values)) {
  test_stat <- tryCatch(ur.df(fred_data_values[[i]], type = "trend", lags = 0)@teststat[1],
    error = function(e) NA
  )
  critical_val <- tryCatch(ur.df(fred_data_values[[i]], type = "trend", lags = 0)@cval,
    error = function(e) NA
  )
  if (!is.na(test_stat) && !is.na(critical_val[2]) && test_stat > critical_val[2]) {
    list_of_non_stationary_trend <- c(list_of_non_stationary_trend, i)
  }
}

cat("Variables non-stationary without trend:\n")
print(unlist(list_of_non_stationary))
cat("Variables non-stationary with trend:\n")
print(unlist(list_of_non_stationary_trend))


# ---------------------- 12. Train-Test Split -------------------
# Define split date: use pre-2020 as training, post-2020 as testing
split_date <- as.Date("2020-01-01")

X_train <- predictors_df %>% filter(sasdate < split_date)
X_test <- predictors_df %>% filter(sasdate >= split_date)

y_train <- cpi_df %>% filter(sasdate < split_date)
y_test <- cpi_df %>% filter(sasdate >= split_date)

cat("Rows: X_train", nrow(X_train), "X_test", nrow(X_test), "y_train", nrow(y_train), "y_test", nrow(y_test), "\n")

summary(X_train$sasdate)
summary(X_test$sasdate)

# ---------------------- 13. Save Cleaned Data -----------------
write_csv(X_train, "../data/fred_train.csv")
write_csv(X_test, "../data/fred_test.csv")
write_csv(y_train, "../data/y_train.csv")
write_csv(y_test, "../data/y_test.csv")



