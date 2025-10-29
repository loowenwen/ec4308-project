library(gbm)

#Data is FX_LEVEL 
data_matrix <- Z_8

parameters_boost = c(5, 10)


normalise <- function(data_matrix) {
  for (col_name in colnames(data_matrix)) {
    col_min <- min(data_matrix[[col_name]], na.rm = TRUE)
    col_max <- max(data_matrix[[col_name]], na.rm = TRUE)
    data_matrix[[col_name]] <- (data_matrix[[col_name]] - col_min) / (col_max - col_min)
  }
  return(data_matrix)
}

boosting_for_one_parameter <- function(train_boosting, val_boosting, parameter){
  mse_matrix = matrix(0, nrow(val_boosting), 1)
  #print("rows_val_boosting")
  #print(nrow(val_boosting))
  for (i in 1:nrow(val_boosting)){
    boostfit = gbm(CPIAUCSL ~., data = train_boosting, distribution='gaussian', bag.fraction=.5,
                   interaction.depth=2, n.trees=parameter, shrinkage=.01)
    #print("HI")
    #print(i)
    #print(parameter)
    value_predicted = predict(boostfit, newdata = val_boosting[i,], n.trees=parameter)
    print(value_predicted)
    mse = (value_predicted - val_boosting$CPIAUCSL[i])^2
    print(mse)
    mse_matrix[i] = mse 
    train_boosting = rbind(train_boosting, val_boosting[i,])
  }
  return(sqrt(mean(mse_matrix, na.rm = TRUE)))
}

gradient_boosting <- function(data_matrix_boost, parameters_boost){
  mse_across_all_parameters = matrix(0, length(parameters_boost))
  boosting_matrix <- data_matrix_boost %>%
    filter(sasdate < as.Date("2020-01-01")) %>%
    select(-sasdate)
  boosting_matrix <- normalise(boosting_matrix)
  train_boosting = boosting_matrix[1:173,]
  val_boosting = boosting_matrix[174:nrow(boosting_matrix),]
  i = 1
  for (j in (parameters_boost)){
    mse_across_all_parameters[i] = boosting_for_one_parameter(train_boosting, val_boosting, j) 
    i = i + 1 
    print(boosting_for_one_parameter(train_boosting, val_boosting, j))
  }
  return(mse_across_all_parameters)
}

print(gradient_boosting(Z_8, parameters_boost))



#mse_across_all_data_all_parameters <- NULL

#for (i in 3:16){
  #data_matrix = get(paste0("Z_",i))
  #mse_across_all_parameters = gradient_boosting(data_matrix, parameters_boost)
  #mse_across_all_data_all_parameters = cbind(mse_across_all_data_all_parameters, mse_across_all_parameters)
#}


