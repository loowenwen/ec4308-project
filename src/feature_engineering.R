library(stats)     
library(zoo)    

# Function to create lagged variables
create_lags <- function(df, lags) {
  lagged <- purrr::map_dfc(1:lags, ~dplyr::lag(df, .x))
  colnames(lagged) <- paste0(rep(names(df), each = lags), "_L", 1:lags)
  return(lagged)
}

##Preparing data to create the various z matrices 
X <- fred_data_clean %>%
  select(-sasdate, -CPIAUCSL) %>%
  select(where(is.numeric)) %>%
  na.omit()



#PCA to create F (take first 25)
pca_result <- prcomp(X, scale. = TRUE)
plot(pca_result, type = "l", main = "Scree Plot of PCA Factors", npcs = 125)
summary(pca_result)
num_factors <- 25
F_t <- as.data.frame(pca_result$x[, 1:num_factors])
colnames(F_t) <- paste0("F", 1:num_factors)

p_y <- 3  #Number of lags for target 
p_f <- 3  #Number of lags for factors 
p_m <- 3 #Number of lags for X

#lagged F (LF_t)
F_lags <- create_lags(F_t, p_f)
#H_t
H_t <- fred_data_values %>% select(-sasdate)
#y_t
y_t <- fred_data_values$CPIAUCSL  
#lagged X_t (LX_t)
X_t_lags <- create_lags(X, p_m)



##F case (LFt)
Z_1 <- bind_cols(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(F_lags))),
  F_lags  
) %>%
  drop_na()

##F-X case
Z_2 <- bind_cols(
  tibble(sasdate = tail(fred_data_values$sasdate, nrow(F_lags))),
  F_lags,
  X_t_lags
) %>% drop_na()

##F-level case
Z_5 <- bind_cols(
  tibble(sasdate = tail(fred_data_values$sasdate, nrow(F_lags))),
  tibble(y_t = tail(y_t, nrow(F_lags))), 
  F_lags,
  tibble(H_t = tail(H_t, nrow(F_lags)))
) %>% drop_na()


##F-X-level case 
Z_8 <- bind_cols(
  tibble(sasdate = tail(fred_data_values$sasdate, nrow(F_lags))),
  F_lags,
  X_t_lags,
  tibble(y_t = tail(y_t, nrow(F_lags))), 
  tibble(H_t = tail(H_t, nrow(F_lags)))
) %>% drop_na()

##X
Z_10 <- bind_cols(
  tibble(sasdate = tail(fred_data_values$sasdate, nrow(X_t_lags))),
  X_t_lags
  ) %>% drop_na()


##X-Level
Z_15 <- bind_cols(
  tibble(sasdate = tail(fred_data_values$sasdate, nrow(X_t_lags))),
  X_t_lags,
  tibble(y_t = tail(y_t, nrow(X_t_lags))), 
  tibble(H_t = tail(H_t, nrow(X_t_lags)))
) %>% drop_na()


##Extra (F, F lag, y, y lag)
y_lag <- create_lags(y_t, p_y)
Z_a <- bind_cols(
  tibble(sasdate = tail(fred_data_values$sasdate, nrow(F_lags))),
  tibble(y_t = tail(y_t, nrow(F_lags))),
  tibble(y_lag = tail(y_lag, nrow(F_lags))),
  tibble(F_t = tail(F_t, nrow(F_lags))),
  F_lags
) %>% drop_na()
