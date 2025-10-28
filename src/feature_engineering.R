library(stats)     
library(zoo)   
library(tidyverse)

# Function to create lagged variables
create_lags <- function(df, p) {
  lags <- lapply(1:p, function(l) dplyr::lag(df, l))
  lagged_df <- do.call(cbind, lags)
  colnames(lagged_df) <- paste0(rep(colnames(df), each = p), "_L", 1:p)
  return(as.data.frame(lagged_df))
}

###MARX function### (general idea: each column is the last p lags of each variable)
create_marx <- function(df, max_lag) {  
  n <- nrow(df)  
  k <- ncol(df)  
  marx_matrix <- matrix(NA, nrow = n - max_lag, ncol = k * max_lag)  
  
  for (p in 1:max_lag) { 
    for (k_idx in 1:k) {  
      col_idx <- (p - 1) * k + k_idx 
      x <- as.matrix(df[, k_idx])
      x_embed <- embed(x, p)
      # Take mean of each row of embedded lags
      x_ma <- rowMeans(x_embed, na.rm = TRUE) / p
      
      # Only keep the last (n - max_lag) rows to align across all p
      marx_matrix[, col_idx] <- tail(x_ma, n - max_lag)
    }
  }
  
  marx_df <- as.data.frame(marx_matrix)
  colnames(marx_df) <- paste0(rep(colnames(df), each = max_lag), "_MA", 1:max_lag)
  
  return(marx_df)
}

###MAF Function###
create_maf <- function(df, n_lags = 12, n_pcs = 1){
  n <- nrow(df)
  k <- ncol(df)
  maf_matrix <- matrix(NA, nrow = n, ncol = k * n_pcs)
  colnames(maf_matrix) <- paste0(rep(colnames(df), each = n_pcs), "_PC", 1:n_pcs)
  
  for (k_idx in 1:k){
    #lagged matrix for this variable
    x <- as.matrix(df[,k_idx])
    lagged_mat <- embed(x, n_lags)   # size of (n - n_lags + 1) x n_lags
    lagged_mat <- na.omit(lagged_mat)
    n_rows <- nrow(lagged_mat)
    
    #PCA
    pca_res <- prcomp(lagged_mat, scale. = TRUE)
    
    #keep first n_pcs components
    pcs <- pca_res$x[, 1:n_pcs, drop = FALSE]
    maf_matrix[(n_lags:n)[1:n_rows], ((k_idx-1)*n_pcs + 1):(k_idx*n_pcs)] <- pcs
  }
  
  return(maf_matrix)
}


###Preparing data to create the various z matrices###
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
#can change lags i just put 3 for now

#lagged F (LF_t)
F_lags <- create_lags(F_t, p_f)

#H_t (non-stationary Xs)
H_t <- fred_data_values %>%
  mutate(sasdate = mdy(sasdate)) %>%  
  filter(year(sasdate) >= 2000 & year(sasdate) <= 2025) %>%  # keep only 2000-2025 rows
  select(-sasdate) %>% 
  select(-CPIAUCSL) %>%
  na.omit()

#y_t
y_t <- fred_data_clean %>%
  na.omit %>% select(-CPIAUCSL)

#X_t
X_t <- X
#lagged X_t (LX_t)
X_t_lags <- create_lags(X, p_m)
#Lagged y_t 
y_lags <- create_lags(y_t, p_y) 

#Feature matrix Z (according to paper - yt, lagged versions of yt, current factors, lagged factors)
Z1_t <- bind_cols(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(X))),
  y_t,
  y_lags,
  F_t,
  F_lags
) %>%
  drop_na()


## Method 2: F-X

#Lag for x
p_x <- 3  # Number of lags for x

#lagged x
x_lags <- create_lags(tibble(x_t = X), p_x)

#Feature matrix Z (according to paper - yt, lagged versions of yt, current factors, lagged factors)
Z2_t <- bind_cols(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(X))),
  y_t,
  y_lags,
  F_t,
  F_lags,
  X_t,
  x_lags
) %>%
  drop_na()

###Setting up the various Z matrices###
#F case (LFt)
Z_1 <- bind_cols(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(F_lags))),
  F_lags  
) %>%
  drop_na()

#F-X case
Z_2 <- bind_cols(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(F_lags))),
  F_lags,
  X_t_lags
) %>% drop_na()


##MARX case (max_lag taken from what was suggested in the paper)
marx_data <- fred_data_clean %>% select(-sasdate) 
Z_11 <- create_marx(marx_data, max_lag = 12) ##MARX with y and x


#F-MARX case
Z_3 <- cbind(
  tail(fred_data_clean$sasdate, nrow(Z_11)),
  tail(F_lags, nrow(Z_11)),
  Z_11
) %>% drop_na()


##MAF Case
p_maf = 12
n_pcs = 2 #just chose 2 but can change acc to optimal oos perf 
maf_data <- fred_data_clean %>% 
  na.omit %>% 
  select(-sasdate) 


Z_12 <- create_maf(maf_data, p_maf, n_pcs) %>% na.omit()

##MAF with y and x

Z_4 <- cbind(
  tail(fred_data_clean$sasdate, nrow(Z_12)),
  tail(F_lags, nrow(Z_12)),
  Z_12
) %>% drop_na()


#F-level case
Z_5 <- bind_cols(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(F_lags))),
  tibble(y_t = tail(y_t, nrow(F_lags))), 
  F_lags,
  tibble(H_t = tail(H_t, nrow(F_lags)))
) %>% drop_na()


#F-X-MARX case
Z_6 <- cbind(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(Z_11))),
  tibble(F_lags = tail(F_lags, nrow(Z_11))),
  tibble(X_t_lags = tail(X_t_lags, nrow(Z_11))),
  Z_11
) %>% drop_na()

#F-X-MAF case
Z_7 <- cbind(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(Z_12))),
  tibble(F_lags = tail(F_lags, nrow(Z_12))),
  tibble(X_t_lags = tail(X_t_lags, nrow(Z_12))),
  Z_12
) %>% drop_na()

##F-X-level case 
Z_8 <- bind_cols(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(F_lags))),
  F_lags,
  X_t_lags,
  tibble(y_t = tail(y_t, nrow(F_lags))), 
  tibble(H_t = tail(H_t, nrow(F_lags)))
) %>% drop_na()

##F-X-MARX-Level
Z_9 <- cbind(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(Z_11))),
  tibble(F_lags = tail(F_lags, nrow(Z_11))),
  tibble(X_t_lags = tail(X_t_lags, nrow(Z_11))),
  Z_11, 
  tibble(H_t = tail(H_t, nrow(Z_11)))
) %>% drop_na()

##X case
Z_10 <- bind_cols(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(X_t_lags))),
  X_t_lags
  ) %>% drop_na()

##X-MARX case
Z_13 <- cbind(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(Z_11))),
  tibble(X_t_lags = tail(X_t_lags, nrow(Z_11))),
  Z_11
) %>% drop_na()

##X-MAF case
Z_14 <- cbind(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(Z_12))),
  tibble(X_t_lags = tail(X_t_lags, nrow(Z_12))),
  Z_12
) %>% drop_na()

##X-Level case
Z_15 <- bind_cols(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(X_t_lags))),
  X_t_lags,
  tibble(y_t = tail(y_t, nrow(X_t_lags))), 
  tibble(H_t = tail(H_t, nrow(X_t_lags)))
) %>% drop_na()

##X-MARX-Level case
Z_16 <- cbind(
  tibble(sasdate = tail(fred_data_clean$sasdate, nrow(Z_11))),
  tibble(X_t_lags = tail(X_t_lags, nrow(Z_11))),
  Z_11,
  tibble(H_t = tail(H_t, nrow(Z_11)))
) %>% drop_na()


