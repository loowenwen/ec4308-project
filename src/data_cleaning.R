###############################################################
# Data Cleaning and Transformation Script
# - Handles missing data
# - Applies FRED-MD transformations (1–7)
# - Forward-fills short gaps (e.g., Apr 2020)
# - Checks for stationarity
# - Splits dataset into training (pre-2020) and testing (post-2020)
###############################################################

# ---------------------- 1. Load Libraries --------------------
library(readr)
library(dplyr)
library(lubridate)
library(zoo)
library(tibble)
library(tidyr)
library(tseries)
library(urca)

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

head(na_summary, 15)  # preview top 15 variables with most missing data

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
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "num_missing") %>%
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
# Extract transformation codes from the first row of the raw file
transform_codes <- as.numeric(fred_raw[1, ])[-1]

# Remove header row to keep actual data
fred_data_values <- fred_clean[-1, ]

# Store transformation functions
transforms <- list(
  transform_1,
  transform_2,
  transform_3,
  transform_4,
  transform_5,
  transform_6,
  transform_7
)

# Apply transformations column by column
for (i in 2:ncol(fred_data_values)) {
  code <- transform_codes[i - 1] # match variable with its code
  fn <- transforms[[code]]       # pick correct transformation
  
  temp <- fn(as.numeric(fred_data_values[[i]]))
  
  # if transformation reduced length, pad with NAs
  if (length(temp) < nrow(fred_data_values)) {
    pad_len <- nrow(fred_data_values) - length(temp)
    temp <- c(rep(NA, pad_len), temp)
  }
  
  fred_data_values[[i]] <- temp
}

# Reapply CPIAUCSL transformation manually (Δlog * 100)
fred_data_values <- fred_data_values %>%
  mutate(CPIAUCSL = 100 * (log(CPIAUCSL) - lag(log(CPIAUCSL))))


# ---------------------- 11. Stationarity Checks ---------------
# Augmented Dickey-Fuller (ADF) test without trend
list_of_non_stationary <- list()

# Loop through each variable (excluding date)
for (i in 2:ncol(fred_data_values)) {
  adf_p <- tryCatch(adf.test(fred_data_values[[i]], k = 0)$p.value,
                    error = function(e) NA)
  if (!is.na(adf_p) && adf_p > 0.05) {
    list_of_non_stationary <- c(list_of_non_stationary, i)
  }
}

# ADF test with trend
list_of_non_stationary_trend <- list()
for (i in 2:ncol(fred_data_values)) {
  test_stat <- tryCatch(ur.df(fred_data_values[[i]], type = "trend", lags = 0)@teststat[1],
                        error = function(e) NA)
  critical_val <- tryCatch(ur.df(fred_data_values[[i]], type = "trend", lags = 0)@cval,
                           error = function(e) NA)
  if (!is.na(test_stat) && !is.na(critical_val[2]) && test_stat > critical_val[2]) {
    list_of_non_stationary_trend <- c(list_of_non_stationary_trend, i)
  }
}

cat("Variables non-stationary without trend:\n")
print(unlist(list_of_non_stationary))
cat("Variables non-stationary with trend:\n")
print(unlist(list_of_non_stationary_trend))


# -- define target variable and macro groups

# inflation target (CPI)
y <- fred_data_values$CPIAUCSL

# ---------------------- 9. Train-Test Split -------------------
# Define split date: use pre-2020 as training, post-2020 as testing
split_date <- as.Date("2020-01-01")

fred_train <- fred_data_values %>%
  filter(sasdate < split_date)

fred_test <- fred_data_values %>%
  filter(sasdate >= split_date)

nrow(fred_train)
nrow(fred_test)

summary(fred_train$sasdate)
summary(fred_test$sasdate)


# ---------------------- 10. Save Cleaned Data -----------------
# Export cleaned and split datasets for feature engineering
write_csv(fred_train, "../data/fred_train.csv")
write_csv(fred_test,  "../data/fred_test.csv")


