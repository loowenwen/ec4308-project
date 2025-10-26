
## RANDOM WALK MODEL ##

fred_data_clean <- read.csv("../data/fred_data_clean.csv")

#Extract Y variable (CPI inflation)
Y <- fred_data_clean$CPIAUCSL

nprev=70 #number of out-of-sample observations (test window )

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


