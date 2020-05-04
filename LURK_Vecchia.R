LURK_Vecchia <- function(Y,X,locs,covparams = NULL,beta.hat = NULL, tol = NULL, m = NULL){
#
#
#
#
#
############################################################## 
### Specify default values for optional inputs  ##############
############################################################## 
N <- length(Y)
if(is.null(covparams)){
  d.sample <- sample(1:N,N/50,replace = FALSE)
  D.sample = rdist(locs[d.sample,1:2])
  t.sample = rdist(locs[d.sample,3])
  covparams <- c(.9*var(Y),mean(D.sample)/4,mean(t.sample)/4,0.1*var(Y)) 
}  
if(is.null(beta.hat)){
  beta.hat <- rep(0,ncol(X))
}
  
if(is.null(tol)){
  tol <- 0.999999
}
  

if(is.null(m)) m = 25

############################################################## 
### Initialize variables ##############
############################################################## 
locs.scaled = cbind(locs[,1]/covparams[2], locs[,2]/covparams[2], locs[,3]/covparams[3])
vecchia.approx=vecchia_specify(locs.scaled,m)
Y.hat <- as.matrix(X)%*%beta.hat


############################################################## 
### Begining algorithm (Algorithm 1 from Messier and Katzfuss 2020) 
############################################################## 
converged=FALSE
prev.error=1e10

while(!converged){
  
  # compute residuals
  res=Y-Y.hat
  vecchia.approx$zord= res[vecchia.approx$ord]
  
  # estimate theta
  vecchia.result<- optim(par=log(covparams),fn=negloglik_vecchia_ST,
                         locs= locs, res=res,vecchia.approx=vecchia.approx,method = "Nelder-Mead",
                         control=list(trace=0))
  
  
  LL.Vecchia.krig <- vecchia.result$value
  covparams=exp(vecchia.result$par)

  
  # transform data to iid
  locs.scaled = cbind(locs[,1]/covparams[2], locs[,2]/covparams[2], locs[,3]/covparams[3])
  vecchia.approx=vecchia_specify(locs.scaled,m)
  
  transformed.data=transform_iid(cbind(Y,as.matrix(X)),
                                 vecchia.approx = vecchia.approx,c(covparams[1],1,0.5),covparams[4])
  
  y.tilde=transformed.data[,1]
  X.tilde=transformed.data[,-1]
  
  # Estimate betas - Vecchia SCAD fitting
  Vecchia.SCAD.fit=cv.ncvreg(as.matrix(X.tilde),y.tilde,family = "gaussian",
                             penalty = "SCAD",dfmax=100,returnX = FALSE)

  idmin <- which(Vecchia.SCAD.fit$lambda == Vecchia.SCAD.fit$lambda.min)
  semin <- Vecchia.SCAD.fit$cve[idmin] + Vecchia.SCAD.fit$cvse[idmin]
  lambda.1se <- max(Vecchia.SCAD.fit$lambda[Vecchia.SCAD.fit$cve<=semin])
  lambda.1se.idx <- which(Vecchia.SCAD.fit$lambda==lambda.1se)
  
  
  ### Betas
  beta.iter <- Vecchia.SCAD.fit$fit$beta[-1,]
  ### Lambdas
  lambda.iter <- Vecchia.SCAD.fit$fit$lambda
  
  ### Get SCAD penalty values
  LL.vecchia.beta <- SCAD_Penalty_Loglike(beta.iter,lambda.iter)
  
  ### Compute log-likelihood
  LL.Vecchia.iter <- LL.Vecchia.krig + LL.vecchia.beta[lambda.1se.idx]
  # Min error (stopping criterion) is the log-likelihood
  min.error=LL.Vecchia.iter

  ### Check min-error against the previous error and tolerance
  if(min.error<prev.error*tol)
  { prev.error=min.error
  beta.hat <- Vecchia.SCAD.fit$fit$beta[,lambda.1se.idx]
  Y.hat <- predict(Vecchia.SCAD.fit,X = X, which = lambda.1se.idx)
  LL.Vecchia <- LL.Vecchia.iter
  covparams.iter <- covparams
  Vecchia.SCAD.iter <- Vecchia.SCAD.fit
  }else{
    converged=TRUE 
    beta.estimates<- beta.hat
    covparams <- covparams.iter
    Vecchia.SCAD.fit <- Vecchia.SCAD.iter
  }
  
}#While (!converged) loop

out.list <- list("covparams" = covparams,"beta" = beta.estimates,"lambda.1se.idx" = lambda.1se.idx, 
                 "Vecchia.SCAD.fit" = Vecchia.SCAD.fit)
return(out.list)

}