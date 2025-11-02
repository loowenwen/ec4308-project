library(dplyr)
library(glmnet)
library(purrr)

---------------------------###USING F-X-MARX (F built using non stationary Xs)####------------------------------------------------
Z_objects <- readRDS("all_Z_matrices.rds")
Za_train <- Z_objects$Z_Fraw_X_MARX_ %>% arrange(sasdate)
Za_test  <- Z_objects$Z_Fraw_X_MARX_test %>% arrange(sasdate)

#Set target: next-period CPI growth
y_train_full <- lead(Za_train$CPI_t_L1)
y_test_full  <- lead(Za_test$CPI_t_L1)

#Drop last row (lead creates NA)
Za_train <- Za_train[-nrow(Za_train), ]
Za_test  <- Za_test[-nrow(Za_test), ]

y_train_full <- y_train_full[-length(y_train_full)]
y_test_full  <- y_test_full[-length(y_test_full)]

#Build the design matrices
feat_train <- Za_train %>% select(-sasdate, -starts_with("CPI_"))
feat_test  <- Za_test  %>% select(-sasdate, -starts_with("CPI_"))

#Keep *only* columns that exist in BOTH sets
common_cols <- intersect(colnames(feat_train), colnames(feat_test))

feat_train <- feat_train %>% select(all_of(common_cols))
feat_test  <- feat_test  %>% select(all_of(common_cols))

x_train_full <- as.matrix(feat_train)
x_test_full  <- as.matrix(feat_test)

#sanity check
stopifnot(ncol(x_train_full) == ncol(x_test_full))
cat("Number of predictors (common):", ncol(x_train_full), "\n")   # should be 1936


MSE <- function(pred, truth) mean((truth - pred)^2, na.rm = TRUE)

# -------------------------------------------------
# Rolling-window Ridge with CV inside each window
# -------------------------------------------------
roll_ridge_cv <- function(x_full, y_full, test_x, test_y,
                          window_size = 120,
                          cv_folds = 10,
                          seed = 1789424) {
  
  n_train <- nrow(x_full)
  n_test  <- nrow(test_x)
  pred_test <- numeric(n_test)
  
  for (t in seq_len(n_test)) {
    # rolling window - expanding to the end of train
    start_idx <- max(1, n_train - window_size + 1)
    end_idx   <- n_train
    
    x_win <- x_full[start_idx:end_idx, , drop = FALSE]
    y_win <- y_full[start_idx:end_idx]
    
    # CV to pick best lambda (only on window data) 
    set.seed(seed)
    cv_fit <- cv.glmnet(x_win, y_win,
                        alpha = 0,          
                        nfolds = cv_folds,
                        standardize = TRUE) 
    
    best_lambda <- cv_fit$lambda.min
    
    #refit on whole window with best lambda 
    ridge_fit <- glmnet(x_win, y_win,
                        alpha = 0,
                        lambda = best_lambda)
    
    # Predict *exactly* the t-th test observation
    pred_test[t] <- predict(ridge_fit,
                            newx = matrix(test_x[t, ], nrow = 1),
                            s = best_lambda)[1, 1]
  }
  
  list(predictions = pred_test,
       mse = MSE(pred_test, test_y))
}

# -------------------------------------------------
# 4. Run the rolling forecast
# -------------------------------------------------
set.seed(1789424)
result <- roll_ridge_cv(x_full = x_train_full,
                        y_full = y_train_full,
                        test_x = x_test_full,
                        test_y = y_test_full,
                        window_size = 120,      # 10 years of monthly data
                        cv_folds = 10)

cat("Rolling-window Ridge MSE on test set:", result$mse, "\n")
head(result$predictions)
