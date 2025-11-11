library(sandwich)
## RANDOM WALK MODEL ##

y_train <- read.csv("../data/y_train.csv")
y_test <- read.csv("../data/y_test.csv")

all_Y <- rbind(y_train, y_test)
#fred_normalized <- data.frame( scale(fred_data_clean[, 2:127]))
#Extract Y variable (CPI inflation)

Y <- all_Y$CPI_t

nprev=100 #number of out-of-sample observations (test window )

oosy=tail(Y,nprev) #get the out-of-sample true values

RMSE <- function(pred, truth){ 
  return(sqrt(mean((truth - pred)^2)))
} 

#create lags
rwtemp=embed(Y,13)
#Simple RW forecast:
rw1c=tail(rwtemp[,2],nprev)       #1-step
rw3c=tail(rwtemp[,4],nprev)       #3-step
rw6c=tail(rwtemp[,7],nprev)       #6-step
rw12c=tail(rwtemp[,13],nprev)     #12-step

#RMSEs for random walk:
rw.rmse1=RMSE(oosy,rw1c)
rw.rmse3=RMSE(oosy,rw3c)
rw.rmse6=RMSE(oosy,rw6c)
rw.rmse12=RMSE(oosy,rw12c)


########################################################
###It's the final showdown: Diebold-Mariano (DM) tests
#########################################################

#####################################################
#Set up for the tests of LASSO with BIC vs. AR benchmark
#####################################################

#Compute squared loss for different horizons (RW)
lrw1c=(oosy-rw1c)^2
lrw3c=(oosy-rw3c)^2
lrw6c=(oosy-rw6c)^2
lrw12c=(oosy-rw12c)^2

# Read in the lasso_bic_results
lasso_bic_results_100_icglm <- readRDS("lasso_bic_results_100_icglm.rds")

# Z_Ht gives the best results for lasso-bic
best_lassobic <- lasso_bic_results_100_icglm$Z_Ht
best_lassobic_1 <- lasso_bic_results_100_icglm$Z_F_X_MARX_Level_naked

#Compute squared loss for different horizons (RF)
lasso1c <- best_lassobic_1$h1$pred
lasso3c <- best_lassobic$h3$pred
lasso6c <- best_lassobic$h6$pred
lasso12c <- best_lassobic$h12$pred

llasso1c=(oosy-lasso1c)^2
llasso3c=(oosy-lasso3c)^2
llasso6c=(oosy-lasso6c)^2
llasso12c=(oosy-lasso12c)^2

#Compute loss differentials (d_t) for different horizons (RW-RF)
drwlasso1=lrw1c-llasso1c
drwlasso3=lrw3c-llasso3c
drwlasso6=lrw6c-llasso6c
drwlasso12=lrw12c-llasso12c

#DM regressions
#Regress d_t (RW-RF) for 1-step forecasts on a constant - get estimate of mean(d_t)
dmrwlasso1=lm(drwlasso1~1) #regression
acf(dmrwlasso1$residuals) #check serial correlation of residuals - number of significant autocorrelations is a good guess for number lags included in the HAC variance estimator
dmrwlasso1$coefficients/sqrt(NeweyWest(dmrwlasso1,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 3-step forecasts on a constant - get estimate of mean(d_t)
dmrwlasso3=lm(drwlasso3~1)
acf(dmrwlasso3$residuals) 
dmrwlasso3$coefficients/sqrt(NeweyWest(dmrwlasso3,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 6-step forecasts on a constant - get estimate of mean(d_t)
dmrwlasso6=lm(drwlasso6~1) 
acf(dmrwlasso6$residuals) 
dmrwlasso6$coefficients/sqrt(NeweyWest(dmrwlasso6,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 12-step forecasts on a constant - get estimate of mean(d_t)
dmrwlasso12=lm(drwlasso12~1) 
acf(dmrwlasso12$residuals) 
dmrwlasso12$coefficients/sqrt(NeweyWest(dmrwlasso12,lag=5)) #form the DM t-statistic

#####################################################
#Set up for the tests of Ridge vs. AR benchmark ###
#####################################################

# Read in the lasso_bic_results
ridge_results_master <- readRDS("ridge_results_100.rds")

# Z_Ht gives the best results for lasso-bic
best_ridge <- ridge_results_master$Z_Ht

#Compute squared loss for different horizons (RF)
ridge1c <- best_ridge$h1$pred
ridge3c <- best_ridge$h3$pred
ridge6c <- best_ridge$h6$pred
ridge12c <- best_ridge$h12$pred

lridge1c=(oosy-ridge1c)^2
lridge3c=(oosy-ridge3c)^2
lridge6c=(oosy-ridge6c)^2
lridge12c=(oosy-ridge12c)^2

#Compute loss differentials (d_t) for different horizons (RW-RF)
drwridge1=lrw1c-lridge1c
drwridge3=lrw3c-lridge3c
drwridge6=lrw6c-lridge6c
drwridge12=lrw12c-lridge12c

#DM regressions
#Regress d_t (RW-RF) for 1-step forecasts on a constant - get estimate of mean(d_t)
dmrwridge1=lm(drwridge1~1) #regression
acf(dmrwridge1$residuals) #check serial correlation of residuals - number of significant autocorrelations is a good guess for number lags included in the HAC variance estimator
dmrwridge1$coefficients/sqrt(NeweyWest(dmrwridge1,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 3-step forecasts on a constant - get estimate of mean(d_t)
dmrwridge3=lm(drwridge3~1)
acf(dmrwridge3$residuals) 
dmrwridge3$coefficients/sqrt(NeweyWest(dmrwridge3,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 6-step forecasts on a constant - get estimate of mean(d_t)
dmrwridge6=lm(drwridge6~1) 
acf(dmrwridge6$residuals) 
dmrwridge6$coefficients/sqrt(NeweyWest(dmrwridge6,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 12-step forecasts on a constant - get estimate of mean(d_t)
dmrwridge12=lm(drwridge12~1) 
acf(dmrwridge12$residuals) 
dmrwridge12$coefficients/sqrt(NeweyWest(dmrwridge12,lag=5)) #form the DM t-statistic

#####################################################
#Set up for the tests of Elstic Net vs. AR benchmark ###
#####################################################
#Compute squared loss for different horizons (EN)
elasticnet_results_100 <- readRDS("elasticnet_results_100.rds")
best_en <- elasticnet_results_100$Z_Ht

en1c <- best_en$h1$pred
en3c <- best_en$h3$pred
en6c <- best_en$h6$pred
en12c <- best_en$h12$pred

len1c=(oosy-en1c)^2
len3c=(oosy-en3c)^2
len6c=(oosy-en6c)^2
len12c=(oosy-en12c)^2

#Compute loss differentials (d_t) for different horizons (ridge-en)
dridgeen1=lridge1c-len1c
dridgeen3=lridge3c-len3c
dridgeen6=lridge6c-len6c
dridgeen12=lridge12c-len12c

#DM regressions
#Regress d_t (RW-RF) for 1-step forecasts on a constant - get estimate of mean(d_t)
dmridgeen1=lm(dridgeen1~1) #regression
acf(dmridgeen1$residuals) #check serial correlation of residuals - number of significant autocorrelations is a good guess for number lags included in the HAC variance estimator
dmridgeen1$coefficients/sqrt(NeweyWest(dmridgeen1,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 3-step forecasts on a constant - get estimate of mean(d_t)
dmrwen3=lm(drwen3~1)
acf(dmrwen3$residuals) 
dmrwen3$coefficients/sqrt(NeweyWest(dmrwen3,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 6-step forecasts on a constant - get estimate of mean(d_t)
dmrwen6=lm(drwen6~1) 
acf(dmrwen6$residuals) 
dmrwen6$coefficients/sqrt(NeweyWest(dmrwen6,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 12-step forecasts on a constant - get estimate of mean(d_t)
dmrwen12=lm(drwen12~1) 
acf(dmrwen12$residuals) 
dmrwen12$coefficients/sqrt(NeweyWest(dmrwen12,lag=5)) #form the DM t-statistic


#####################################################
#Set up for the tests of Ensemble with AR benchmark
#####################################################

# Read in the ensemble results 
ensemble_best <- readRDS("all_best_models_preds.rds")
ensemble_grunrc <- ensemble_best$grunrc_ensemble

#Compute squared loss for different horizons (GRUNRC)
grunrc1c <- ensemble_grunrc$h1$pred
grunrc3c <- ensemble_grunrc$h3$pred
grunrc6c <- ensemble_grunrc$h6$pred
grunrc12c <- ensemble_grunrc$h12$pred

lgrunrc1c=(oosy-grunrc1c)^2
lgrunrc3c=(oosy-grunrc3c)^2
lgrunrc6c=(oosy-grunrc6c)^2
lgrunrc12c=(oosy-grunrc12c)^2

#Compute loss differentials (d_t) for different horizons (RW-GRUNRC)
drwgrunrc1=lrw1c-lgrunrc1c
drwgrunrc3=lrw3c-lgrunrc3c
drwgrunrc6=lrw6c-lgrunrc6c
drwgrunrc12=lrw12c-lgrunrc12c

#DM regressions
#Regress d_t (RW-GRUNRC) for 1-step forecasts on a constant - get estimate of mean(d_t)
dmrwgrunrc1=lm(drwgrunrc1~1) #regression
acf(dmrwgrunrc1$residuals) #check serial correlation of residuals - number of significant autocorrelations is a good guess for number lags included in the HAC variance estimator
dmrwgrunrc1$coefficients/sqrt(NeweyWest(dmrwgrunrc1,lag=5)) #form the DM t-statistic

#Regress d_t (RW-GRUNRC) for 1-step forecasts on a constant - get estimate of mean(d_t)
dmrwgrunrc3=lm(drwgrunrc3~1) #regression
acf(dmrwgrunrc3$residuals) #check serial correlation of residuals - number of significant autocorrelations is a good guess for number lags included in the HAC variance estimator
dmrwgrunrc3$coefficients/sqrt(NeweyWest(dmrwgrunrc3,lag=5)) #form the DM t-statistic

#Regress d_t (RW-GRUNRC) for 6-step forecasts on a constant - get estimate of mean(d_t)
dmrwgrunrc6=lm(drwgrunrc6~1) #regression
acf(dmrwgrunrc6$residuals) #check serial correlation of residuals - number of significant autocorrelations is a good guess for number lags included in the HAC variance estimator
dmrwgrunrc6$coefficients/sqrt(NeweyWest(dmrwgrunrc6,lag=5)) #form the DM t-statistic

#Regress d_t (RW-GRUNRC) for 12-step forecasts on a constant - get estimate of mean(d_t)
dmrwgrunrc12=lm(drwgrunrc12~1) #regression
acf(dmrwgrunrc12$residuals) #check serial correlation of residuals - number of significant autocorrelations is a good guess for number lags included in the HAC variance estimator
dmrwgrunrc12$coefficients/sqrt(NeweyWest(dmrwgrunrc12,lag=5)) #form the DM t-statistic

#DM regressions - Comparing Best 2 Models for each Horizon 
#Regress d_t (GRUNRC-Ridge) for 1-step forecasts on a constant - get estimate of mean(d_t)
dgrunrclasso1 = lgrunrc1c - lridge1c
dmgrunrclasso1=lm(dgrunrclasso1~1) 
acf(dmgrunrclasso1$residuals) 
dmgrunrclasso1$coefficients/sqrt(NeweyWest(dmgrunrclasso1,lag=5)) #form the DM t-statistic

#Regress d_t (LASSO-Ridge) for 3-step forecasts on a constant - get estimate of mean(d_t)
dlassoridge3=lridge3c-llasso3c
dmlassoridge3=lm(dlassoridge3~1)
acf(dmlassoridge3$residuals) 
dmlassoridge3$coefficients/sqrt(NeweyWest(dmlassoridge3,lag=5)) #form the DM t-statistic

#Regress d_t (GRUNRC-Ridge) for 6-step forecasts on a constant - get estimate of mean(d_t)
dgrunrclasso6 = lgrunrc6c - lridge6c
dmgrunrclasso6=lm(dgrunrclasso6~1) 
acf(dmgrunrclasso6$residuals) 
dmgrunrclasso6$coefficients/sqrt(NeweyWest(dmgrunrclasso6,lag=5)) #form the DM t-statistic

#Regress d_t (RF-EN) for 12-step forecasts on a constant - get estimate of mean(d_t)
denlasso12=len12c-llasso12c
dmenlasso12=lm(denlasso12~1) 
acf(dmenlasso12$residuals) 
dmenlasso12$coefficients/sqrt(NeweyWest(dmenlasso12,lag=5)) #form the DM t-statistic


#Regress d_t (RW-RF) for 1-step forecasts on a constant - get estimate of mean(d_t)
# dmrflasso6=lm(drflasso6~1) #regression
# acf(dmrflasso6$residuals) #check serial correlation of residuals - number of significant autocorrelations is a good guess for number lags included in the HAC variance estimator
# dmrflasso6$coefficients/sqrt(NeweyWest(dmrflasso6,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 1-step forecasts on a constant - get estimate of mean(d_t)
#dmenlasso12=lm(denlasso12~1) #regression
#acf(dmenlasso12$residuals) #check serial correlation of residuals - number of significant autocorrelations is a good guess for number lags included in the HAC variance estimator
#dmenlasso12$coefficients/sqrt(NeweyWest(dmenlasso12,lag=5)) #form the DM t-statistic

#Regress d_t (RW-RF) for 6-step forecasts on a constant - get estimate of mean(d_t)
#dmrfen6=lm(drfen6~1) 
#acf(dmrfen6$residuals) 
#dmrfen6$coefficients/sqrt(NeweyWest(dmrfen6,lag=5)) #form the DM t-statistic
#Regress d_t (GRUNRC-Ridge) for 3-step forecasts on a constant - get estimate of mean(d_t)
# dmrfen3=lm(drfen3~1)
# acf(dmrfen3$residuals) 
# dmrfen3$coefficients/sqrt(NeweyWest(dmrfen3,lag=5)) #form the DM t-statistic

#####################################################################################################
#Compare Random Forest vs LASSO
#####################################
# random_forest_master <- readRDS("../data/rf_results/rf_master_results.rds")
# 
# # Z_Ht gives the best results for lasso-bic
# best_rf <- ridge_results_master$Z_Ht
# 
# #Compute squared loss for different horizons (RF)
# rf1c <- best_rf$h1$pred
# rf3c <- best_rf$h3$pred
# rf6c <- best_rf$h6$pred
# rf12c <- best_rf$h12$pred
# 
# lrf1c=(oosy-rf1c)^2
# lrf3c=(oosy-rf3c)^2
# lrf6c=(oosy-rf6c)^2
# lrf12c=(oosy-rf12c)^2
# 
# drflasso6=lrf6c-llasso6c

#####################################################################################################
#Compare Ridge vs LASSO
#####################################
#Compute loss differentials (d_t) for different horizons (RW-RF)
# dlassoridge1=lridge1c-llasso1c
# dlassoridge3=lridge3c-llasso3c
# dlassoridge6=lridge6c-llasso6c
# dlassoridge12=lridge12c-llasso12c

