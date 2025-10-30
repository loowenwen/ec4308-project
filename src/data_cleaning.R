## This script prepares the FRED-MD macroeconomic dataset for analysis by 
## applying the official transformation codes, cleaning the data, checking for 
## stationarity, and organizing variables by category.

# -- setup and library loading --

# load required libraries
library(readr)
library(dplyr)
library(lubridate)
library(tseries)
library(urca)
library(zoo)

# -- load raw FRED-MD data --

# read dataset
fred_data <- read_csv("../data/2025-09-MD.csv")
head(fred_data)

## The first row contains transformation codes, while subsequent rows contain 
## the actual time series values.


# -- define transformation functions --

# code 1: no transformation (level)
transform_1 <- function(x) {
  x
}

# code 2: first difference
transform_2 <- function(x) {
  x - lag(x, n = 1, default = NA)
}

# code 3: second difference
transform_3 <- function(x) {
  (x - lag(x, 1)) - (lag(x, 1) - lag(x, 2))
}

# code 4: log transform
transform_4 <- function(x) {
  ifelse(x > 0, log(x), NA) # safer than plain log(x)
}

# code 5: log first difference
transform_5 <- function(x) {
  log(x) - lag(log(x), 1)
}

# code 6: log second difference
transform_6 <- function(x) {
  (log(x) - lag(log(x), 1)) - (lag(log(x), 1) - lag(log(x), 2))
}

# code 7: custom ratio transform (rarely used)
transform_7 <- function(x) {
  (x / (lag(x, 1) - 1)) - (lag(x, 1) / (lag(x, 2) - 1))
}


# -- apply transformations according to FRED-MD codes

# extract transformation codes from the first row (excluding date)
transform_codes <- as.numeric(fred_data[1, ])[-1]

# store actual data (from row 2 onwards)
fred_data_values <- fred_data[-1, ]

# create list of transformation functions
transforms <- list(
  transform_1,
  transform_2,
  transform_3,
  transform_4,
  transform_5,
  transform_6,
  transform_7
)

# apply transformations column by column
# pad with NAs if differencing shortens the vector
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


# -- check number of NAs per column --

# count NAs in each column
na_count <- fred_data_values %>%
  summarise(across(everything(), ~ sum(is.na(.))))

# view the result in an easy-to-read format
na_count_long <- na_count %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "num_missing") %>%
  arrange(desc(num_missing))

head(na_count_long, 15)

# -- check correlation of variables that have high number of NAs --
## namely "ACOGNO", "UMCSENTx", "TWEXAFEGSMTHx", "ANDENOx"

# correlation testing of "ACOGNO"
target_var <- "ACOGNO"

# compute correlation with all other numeric columns
corrs <- fred_data_values %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(),
                   ~ cor(.x, fred_data_values[[target_var]], use = "pairwise.complete.obs"))) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "correlation") %>%
  arrange(desc(abs(correlation)))

head(corrs, 10)

# correlation testing of "UMCSENTx"
target_var <- "UMCSENTx"

# compute correlation with all other numeric columns
corrs <- fred_data_values %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(),
                   ~ cor(.x, fred_data_values[[target_var]], use = "pairwise.complete.obs"))) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "correlation") %>%
  arrange(desc(abs(correlation)))

head(corrs, 10)

# correlation testing of "TWEXAFEGSMTHx"
target_var <- "TWEXAFEGSMTHx"

# compute correlation with all other numeric columns
corrs <- fred_data_values %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(),
                   ~ cor(.x, fred_data_values[[target_var]], use = "pairwise.complete.obs"))) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "correlation") %>%
  arrange(desc(abs(correlation)))

head(corrs, 10)

# correlation testing of "ANDENOx"
target_var <- "ANDENOx"

# compute correlation with all other numeric columns
corrs <- fred_data_values %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(),
                   ~ cor(.x, fred_data_values[[target_var]], use = "pairwise.complete.obs"))) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "correlation") %>%
  arrange(desc(abs(correlation)))

head(corrs, 10)

# -- remove redundant columns --

fred_data_clean <- fred_data_values %>%
  select(-ACOGNO) %>%        # remove ACOGNO
  mutate(sasdate = mdy(sasdate)) %>%  # ensure date is in Date format
  drop_na(UMCSENTx)           # drop only rows where UMCSENTx is NA

na_summary <- fred_data_clean %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "num_missing") %>%
  arrange(desc(num_missing))

print(na_summary)

fred_data_clean %>%
  filter(if_any(everything(), is.na))

# calculate monthly gaps in days
fred_data_clean %>%
  arrange(sasdate) %>%
  mutate(diff_days = as.numeric(sasdate - lag(sasdate))) %>%
  summarise(
    total_obs = n(),
    min_gap = min(diff_days, na.rm = TRUE),
    max_gap = max(diff_days, na.rm = TRUE),
    missing_gaps = sum(diff_days > 40, na.rm = TRUE)   # more than ~1 month
  )

fred_data_clean %>%
  arrange(sasdate) %>%
  mutate(diff_days = as.numeric(sasdate - lag(sasdate))) %>%
  filter(diff_days > 40)

## 1 2020-04-01 2 2020-05-01 3 2025-07-01 4 2025-08-01

# define the dates of interest
target_dates <- as.Date(c("2020-04-01", "2020-05-01", "2025-07-01", "2025-08-01"))

# filter only those rows
missing_rows <- fred_data_clean %>%
  filter(sasdate %in% target_dates)

# for each of these rows, find which columns are NA
missing_map <- missing_rows %>%
  pivot_longer(cols = -sasdate,
               names_to = "variable",
               values_to = "value") %>%
  filter(is.na(value)) %>%
  arrange(sasdate)

print(missing_map)

# -- handle the missing NAs in these 4 rows --

# remove incomplete future months 
fred_data_clean <- fred_data_clean %>%
  filter(sasdate <= as.Date("2025-06-01"))   # keep up to June 2025 only

# impute missing values in test period only
fred_data_clean <- fred_data_clean %>%
  arrange(sasdate) %>%
  mutate(across(where(is.numeric),
                ~ na.approx(.x, na.rm = FALSE)))   # linear interpolation

# check 
fred_data_clean %>%
  filter(sasdate %in% target_dates) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "num_missing") %>%
  filter(num_missing > 0)


# -- stationarity check with Augmented Dickey-Fuller (ADF) tests
list_of_non_stationary <- list()

# Loop through each variable (excluding date)
for (i in 2:ncol(fred_data_clean)) {
  adf_p <- tryCatch(adf.test(fred_data_clean[[i]], k = 0)$p.value,
                    error = function(e) NA)
  if (!is.na(adf_p) && adf_p > 0.05) {
    list_of_non_stationary <- c(list_of_non_stationary, i)
  }
}

# ADF with trend
list_of_non_stationary_trend <- list()
for (i in 2:ncol(fred_data_clean)) {
  test_stat <- tryCatch(ur.df(fred_data_clean[[i]], type = "trend", lags = 0)@teststat[1],
                        error = function(e) NA)
  critical_val <- tryCatch(ur.df(fred_data_clean[[i]], type = "trend", lags = 0)@cval,
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
y <- fred_data_clean$CPIAUCSL

which(names(fred_data) == "ACOGNO")

# example predictor groupings based on FRED-MD classification
x_cpi <- fred_data_clean[, 106:108]
x_g1_outputincome <- fred_data_clean[, c(2:3, 7:20)]
x_g2_labormkt <- fred_data_clean[, c(21:48, 119:121)]
x_g3_consumptionorders <- fred_data_clean[, 49:58]
x_g4_ordersinv <- fred_data_clean[, c(4:6, 59:63, 122)]
x_g5_moneycredit <- fred_data_clean[, c(64:73, 124:125)]
x_g6_intexrates <- fred_data_clean[, 77:98]
x_g7_prices <- fred_data_clean[, 103:118]
x_g8_stockmkt <- fred_data_clean[, 74:76]

# identify unused variables (optional)
all_used <- c(107:109, 2:3, 7:20, 21:48, 120:122,
              49:58, 4:6, 59:64, 123,
              65:74, 125:126, 78:99,
              104:119, 75:77)
x_other <- fred_data_clean[, -all_used]

# -- save clean dataset --
train_data <- fred_data_clean %>%
  filter(sasdate <= as.Date("2019-12-01"))

test_data <- fred_data_clean %>%
  filter(sasdate >= as.Date("2020-01-01"))

nrow(train_data)
nrow(test_data)

summary(train_data$sasdate)
summary(test_data$sasdate)

# save training set
write_csv(train_data, "../data/fred_train.csv")

# save test set
write_csv(test_data, "../data/fred_test.csv")
