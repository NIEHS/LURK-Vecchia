
###### Libraries needed in Vecchia ###########
library(GPvecchia)
library(Matrix)
library(fields)
### SCAD
library(ncvreg)

# For reading in Excel data 
library(readxl)
# Proper scoring rules (log-score and CRPS)
library(scoringRules)


source("U:/Papers/Stats_Paper/varSelection_kyle/GPvecchia_Model_Selection/lurk-vecchia/LURK_Vecchia.R")
source("U:/Papers/Stats_Paper/varSelection_kyle/GPvecchia_Model_Selection/lurk-vecchia/LURK_Helpers.R")

###### Working directory
setwd("U:/Papers/Stats_Paper/varSelection_kyle/US_NO2/")


### DATA #####
US_NO2_Data<- read.csv("US_NO2_uniqueST_wCovariates_20191001.csv")

lat       <- US_NO2_Data$Latitude # Latitude
lon       <- US_NO2_Data$Longitude # Longitude
logNO2     <-  log(US_NO2_Data$Y) # ozone data
time      <- US_NO2_Data$YearFrac    # time in fractional years
N         <- length(logNO2)
locs = cbind(lon,lat,time) # Coordinates in R3 (x,y,t)

X.Design <- US_NO2_Data[,7:145] #Design matrix
X.Design <- scale(X.Design)
ID <- US_NO2_Data$ID
# Set the seed and create foldids - seed is necessary for repeatability
set.seed(1)
foldid <- sample(1:10,N,replace=TRUE)

Local.Kriging.MSE <- matrix(nrow = 10)
Local.Kriging.crps <- matrix(nrow = 10)
Local.Kriging.logs <- matrix(nrow = 10)

LUR.iid.MSE <- matrix(nrow = 10)
LUR.iid.crps <- matrix(nrow = 10)
LUR.iid.logs <- matrix(nrow = 10)

LURK.Local.MSE <- matrix(nrow = 10)
LURK.Local.crps <- matrix(nrow = 10)
LURK.Local.logs <- matrix(nrow = 10)


LURK.Vecchia.MSE <- matrix(nrow = 10)
LURK.Vecchia.crps <- matrix(nrow = 10)
LURK.Vecchia.logs <- matrix(nrow = 10)


########################
### Begin Main Loop - 10 fold cross-validation of NO2 dataset #####
for (i in 1:10){
  
  
  # Set up the training and test data
  idx.train <- foldid != i
  idx.test  <- foldid == i
  
  #### FOR TESTING PURPOSE - use test set size
  # idx.train <- idx.test
  ###
  locs.train <- locs[idx.train,]
  locs.test  <- locs[idx.test,]
  X.train <- X.Design[idx.train,]
  X.test  <- X.Design[idx.test,]
  Y.train <- logNO2[idx.train]
  Y.test  <- logNO2[idx.test]
  n=sum(idx.train)
  
  ################### 
  ### Pure Geostatistical estimation (Local Simple Kriging)
  ###################
  print(i)
  print("OK")
  
  
  OK.cov.param = ST_Krig_Param_Avg(Y.train,locs.train,p = 500,k = 10)
  
  #(new_coords,obs_coords,Y_obs,cov.pars,NN)
  OK.pred <- Kr_pred(locs.test,locs.train,Y.train,OK.cov.param,25)  
  OK.mean<- OK.pred$Kr.prediction
  OK.sd <- sqrt(OK.pred$Kr.Var + OK.cov.param[4])
  
  Local.Kriging.MSE[i] <- mean((OK.mean- Y.test)^2)
  Local.Kriging.logs[i] <- mean(logs_norm(Y.test,OK.mean,OK.sd))    
  Local.Kriging.crps[i] <- mean(crps_norm(Y.test,OK.mean,OK.sd))    
  
  ######## iid SCAD #############
  print(i)
  print("SCAD")
  SCAD.fit=cv.ncvreg(as.matrix(X.train),Y.train,family = "gaussian",
                     penalty = "SCAD",returnX = FALSE)
  
  idmin <- which(SCAD.fit$lambda == SCAD.fit$lambda.min)
  semin <- SCAD.fit$cve[idmin] + SCAD.fit$cvse[idmin]
  lambda.1se <- max(SCAD.fit$lambda[SCAD.fit$cve<=semin])
  lambda.1se.idx <- which(SCAD.fit$lambda==lambda.1se)
  
  SCAD.pred <- predict(SCAD.fit,X = as.matrix(X.test),lambda=lambda.1se) 
  SCAD.beta.estimates <- coef(SCAD.fit,lambda=lambda.1se)
  beta.in <- SCAD.beta.estimates!=0
  # Get the SCAD prediction variance assuming LR assumptions  
  tau2.SCAD <- (1/(n-sum(SCAD.beta.estimates!=0))) * norm(Y.test - SCAD.pred,type = "F")^2
  
  LUR.iid.MSE[i] = mean((SCAD.pred-Y.test)^2)
  LUR.iid.crps[i] <- mean(crps_norm(Y.test,SCAD.pred, rep(sqrt(tau2.SCAD),length(SCAD.pred))))
  LUR.iid.logs[i] <- mean(logs_norm(Y.test,SCAD.pred, rep(sqrt(tau2.SCAD),length(SCAD.pred))))
  
  ######## Kriging with an external drift  #############
  # de-trend
  print(i)
  print("KED")
  SCAD.train.pred <- predict(SCAD.fit,X = as.matrix(X.train),lambda=lambda.1se)
  res.train <- Y.train- SCAD.train.pred
  
  KED.cov.param = ST_Krig_Param_Avg(res.train,locs.train,p = 500,k = 10)
  
  
  KED.pred <- Kr_pred(locs.test,locs.train,res.train,KED.cov.param,25)  
  KED.mean<- KED.pred$Kr.prediction + SCAD.pred
  KED.sd <- sqrt(KED.pred$Kr.Var + KED.cov.param[4])
  
  
  
  LURK.Local.MSE[i] <- mean((KED.mean- Y.test)^2)
  LURK.Local.logs[i] <- mean(logs_norm(Y.test,KED.mean,KED.sd))    
  LURK.Local.crps[i] <- mean(crps_norm(Y.test,KED.mean,KED.sd))    
  
  
  
  ##############################################
  #########LURK-Vecchia ########################
  ##############################################
  
  m=25
  results.LURK.Vecchia <- LURK_Vecchia(Y = Y.train,X = X.train,locs = locs.train, m = m)
  
  
 
  Vecchia.beta.estimates <- results.LURK.Vecchia$beta  
  Vecchia.theta <- results.LURK.Vecchia$covparam
  
  # Vecchia prediction at new locations
  LURK.Vecchia.Pred <- predict(results.LURK.Vecchia$Vecchia.SCAD.fit,X = X.test,which = results.LURK.Vecchia$lambda.1se.idx) 
  # Vecchia trend prediction at observed data
  LURK.Vecchia.hat <- predict(results.LURK.Vecchia$Vecchia.SCAD.fit,X = X.train,which = results.LURK.Vecchia$lambda.1se.idx)
  
  # Test set prediction
  res=Y.train- LURK.Vecchia.hat
  
  locs.train.scaled = cbind(locs.train[,1]/Vecchia.theta[2], 
        locs.train[,2]/Vecchia.theta[2], locs.train[,3]/Vecchia.theta[3])  
  
  locs.test.scaled = cbind(locs.test[,1]/Vecchia.theta[2], 
                           locs.test[,2]/Vecchia.theta[2], locs.test[,3]/Vecchia.theta[3])  
  
  
  vec.approx.test=vecchia_specify(locs.train.scaled,m,locs.pred=locs.test.scaled)
  

  ## carry out prediction
  pred=vecchia_prediction(res,vec.approx.test,c(Vecchia.theta[1],1,0.5),Vecchia.theta[4])
  
  # returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord
  Vec.mean=pred$mu.pred+LURK.Vecchia.Pred
  Vec.sds=sqrt(pred$var.pred + Vecchia.theta[4])
  
  LURK.Vecchia.MSE[i] <- mean((Vec.mean - Y.test)^2)
  LURK.Vecchia.logs[i] <- mean(logs_norm(Y.test,Vec.mean,Vec.sds))    
  LURK.Vecchia.crps[i] <- mean(crps_norm(Y.test,Vec.mean,Vec.sds))   
  
  
  
}#Outer loop (10 fold cross-validation)



save(Local.Kriging.crps,Local.Kriging.logs,Local.Kriging.MSE,
     LUR.iid.crps,LUR.iid.logs,LUR.iid.MSE,
     LURK.Local.crps,LURK.Local.logs,LURK.Local.MSE,
  LURK.Vecchia.crps,LURK.Vecchia.logs,LURK.Vecchia.MSE,
      file = "ValStats_NO2_10fold_Crossvalidation.RData")


