library(stats)     
library(zoo)    

#y and x matrices
y <- fred_data_clean$CPIAUCSL
X <- fred_data_clean %>%
  select(-sasdate, -CPIAUCSL) %>%
  select(where(is.numeric)) %>%
  na.omit()

# Match y to X row count after cleaning
y <- tail(y, nrow(X))

#PCA
pca_result <- prcomp(X, scale. = TRUE)
plot(pca_result, type = "l", main = "Scree Plot of PCA Factors", npcs = 125)
summary(pca_result)
num_factors <- 25

# first 25 PCs
F_t <- as.data.frame(pca_result$x[, 1:num_factors])
colnames(F_t) <- paste0("F", 1:num_factors)

#Lags for y and factors
p_y <- 3  # Number of lags for target 
p_f <- 3  # Number of lags for factors 

# Function to create lagged variables
create_lags <- function(df, lags) {
  lagged <- purrr::map_dfc(1:lags, ~dplyr::lag(df, .x))
  colnames(lagged) <- paste0(rep(names(df), each = lags), "_L", 1:lags)
  return(lagged)
}

#lagged y
y_lags <- create_lags(tibble(y_t = y), p_y)

#lagged factors (PCs)
F_lags <- create_lags(F_t, p_f)

#Feature matrix Z (according to paper - yt, lagged versions of yt, current factors, lagged factors)
Z_t <- bind_cols(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(X))),
  tibble(y_t = y),
  y_lags,
  F_t,
  F_lags
) %>%
  drop_na()

head(Z_t)  
dim(Z_t)  