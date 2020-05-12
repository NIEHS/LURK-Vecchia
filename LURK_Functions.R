################################################################################
### LURK-Vecchia: Spatiotemporal Land-use Regression Kriging with Vecchia Approximation
################################################################################

LURK_Vecchia <- function(Y,X,locs,covparams = NULL,beta.hat = NULL, tol = NULL, m = NULL){

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


###########################################################
## LURK-Full ##############################################
###############################################

LURK_Full <- function(Y,X,locs,covparams = NULL,beta.hat = NULL, tol = NULL){
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
  Y.hat <- as.matrix(X)%*%beta.hat
  
  
  ############################################################## 
  ### Begining algorithm (Algorithm 1 from Messier and Katzfuss 2020, 
  #      but for the Full model) 
  ############################################################## 
  converged=FALSE
  prev.error=1e10
  
  while(!converged){
    
    # compute residuals
    res=Y-Y.hat
    
    # estimate theta
    
    full.result=optim(par=log(covparams),fn=negloglik_full_ST,
                      locs= locs, res=res,N=N,method = "Nelder-Mead",
                      control=list(trace=0))
    
    
    covparams <- exp(full.result$par)
    LL.krig <- full.result$value
    
    # transform data to iid
    locs.scaled = cbind(locs[,1]/covparams[2], locs[,2]/covparams[2], locs[,3]/covparams[3])
    
    Omega.full <- full.theta.hat[1]*Exponential(rdist(locs.scaled),
                                                range=1)+full.theta.hat[4]*diag(N)
    Omega.lc <- solve(t(chol(Omega.full)))
    
    
    y.tilde=transformed.data[,1]
    X.tilde=transformed.data[,-1]
    
    # Estimate betas - Full SCAD fitting
    Full.SCAD.fit=cv.ncvreg(as.matrix(X.tilde),y.tilde,family = "gaussian",
                            penalty = "SCAD",dfmax=100,returnX = FALSE)
    
    idmin <- which(full.SCAD.fit$lambda == full.SCAD.fit$lambda.min)
    semin <- full.SCAD.fit$cve[idmin] + full.SCAD.fit$cvse[idmin]
    lambda.1se <- max(full.SCAD.fit$lambda[full.SCAD.fit$cve<=semin])
    lambda.1se.idx <- which(full.SCAD.fit$lambda==lambda.1se)
    
    
    ### Betas
    beta.iter <- Full.SCAD.fit$fit$beta[-1,]
    ### Lambdas
    lambda.iter <- Full.SCAD.fit$fit$lambda
    
    ### Get SCAD penalty values
    LL.full.beta <- SCAD_Penalty_Loglike(beta.iter,lambda.iter)
    
    ### Compute log-likelihood
    LL.full.iter <- LL.krig + LL.full.beta[lambda.1se.idx]
    # Min error (stopping criterion) is the log-likelihood
    min.error=LL.full.iter
    
    ### Check min-error against the previous error and tolerance
    if(min.error<prev.error*tol)
    { prev.error=min.error
    beta.hat <- Full.SCAD.fit$fit$beta[,lambda.1se.idx]
    Y.hat <- predict(Full.SCAD.fit,X = X, which = lambda.1se.idx)
    LL.Full <- LL.full.iter
    covparams.iter <- covparams
    Full.SCAD.iter <- Full.SCAD.fit
    }else{
      converged=TRUE 
      beta.estimates<- beta.hat
      covparams <- covparams.iter
      Full.SCAD.fit <- Full.SCAD.iter
    }
    
  }#While (!converged) loop
  
  out.list <- list("covparams" = covparams,"beta" = beta.estimates,"lambda.1se.idx" = lambda.1se.idx, 
                   "Full.SCAD.fit" = Full.SCAD.fit)
  return(out.list)
  
}


################################################################################
# Spatiotemporal Vecchia negative loglikelihood
################################################################################

negloglik_vecchia_ST=function(logparms,locs,res,vecchia.approx){
  parms = exp(logparms)
  locs.scaled = cbind(locs[,1]/parms[2], locs[,2]/parms[2], locs[,3]/parms[3])
  vecchia.approx$locsord = locs.scaled[vecchia.approx$ord,]
  -vecchia_likelihood(res, vecchia.approx, c(parms[1],1, 0.5), parms[4] ) 
}  


################################################################################
# Spatiotemporal Full Kriging negative loglikelihood
################################################################################
negloglik_full_ST=function(logparms,locs,y,N){
  parms = exp(logparms)
  locs.scaled = cbind(locs[,1]/parms[2], locs[,2]/parms[2], locs[,3]/parms[3])
  d <- fields::rdist(locs.scaled)
  cov.mat=parms[1]*fields::Exponential(d,range=1)+
    parms[4]*diag(N)
  -mvtnorm::dmvnorm(y,rep(0,N),cov.mat,log=TRUE)
}


################################################################################
# SCAD Penalty value
################################################################################
SCAD_Penalty_Loglike <- function(beta.in,lambda){
  
  penalty <- matrix(NA,nrow = nrow(beta.in),ncol = ncol(beta.in))
  for (j in 1:nrow(beta.in)){
    beta.j <- beta.in[j,]
    idx1 <- abs(beta.j)<= lambda
    idx2 <- abs(beta.j)>lambda & abs(beta.j)<=3.7*lambda
    idx3 <- abs(beta.j)>3.7*lambda
    penalty[j,idx1] <- lambda[idx1]*beta.j[idx1]
    penalty[j,idx2] <- -(abs(beta.j[idx2])^2 - 7.4*lambda[idx2]*abs(beta.j[idx2])+lambda[idx2]^2)/(5.4)
    penalty[j,idx3] <- (3.7*lambda[idx3]^2 + lambda[idx3]^2)/2
  }
  
  
  loglik.penalty <- lambda*colSums(penalty)
  
  return(loglik.penalty)
}
# 
################################################################################
# Ordinary Spatiotemporal Kriging Prediction with a Local S-T neigborhood
################################################################################
Kr_pred  =  function(new_coords,obs_coords,Y_obs,cov.pars,NN){
  
  Kr.prediction = matrix(NA,nrow = nrow(new_coords),ncol=1)
  Kr.Var = matrix(NA,nrow = nrow(new_coords),ncol=1)
  
  new_coords_scaled <- cbind(new_coords[,1]/cov.pars[2],new_coords[,2]/cov.pars[2],new_coords[,3]/cov.pars[3])
  obs_coords_scaled <- cbind(obs_coords[,1]/cov.pars[2],obs_coords[,2]/cov.pars[2],obs_coords[,3]/cov.pars[3])
  df_new <- new_coords_scaled
  df_obs <- obs_coords_scaled
  
  for (i in 1:nrow(new_coords)){
   # print(i)
    locs.test.i <- rep.row(df_new[i,],nrow(df_obs))
    dist.op<- fields::rdist.vec(locs.test.i,df_obs)
    Sigma.op <- cov.pars[1] * fields::Exponential(dist.op,range = 1)
    s.op <- sort(dist.op,index.return = TRUE)
    s.idx <- s.op$ix[1:NN]
    # Get the closest "NN" observations
    dist.oo <- fields::rdist(df_obs[s.idx,],df_obs[s.idx,])
    Sigma.oo <- cov.pars[4]*diag(NN) + 
      cov.pars[1]*fields::Exponential(dist.oo,range=1)
    oo <- Y_obs[s.idx]
    ## kriging predictor (posterior mean)
    mean_trend <- mean(oo)
    Kr.prediction[i] <- mean_trend + t(Sigma.op[s.idx]) %*% solve(Sigma.oo,oo - mean_trend)
    
    Kr.Var[i] = cov.pars[1]-t(Sigma.op[s.idx]) %*% solve(Sigma.oo, Sigma.op[s.idx])
  }
  
  Kr.Data <- data.frame(Kr.prediction,Kr.Var)
  
  return(Kr.Data)
}


################################################################################
# Spatiotemporal Kriging Maximum Likelihood Estiamtion based on an average
# of multiple random subsets
################################################################################

ST_Krig_Param_Avg = function(Y,locs,p,k = 10){
  

n = length(Y)
mdl.geo.fit.avg <- list()
mdl.geo.fit.avg <- matrix(NA,nrow=k,ncol = 4)

for (i in 1:k){
  set.seed(i)
  OK.fold = sample(1:n,p,replace=FALSE)
  Y.train = Y[OK.fold]
  locs.train <- locs[OK.fold,]
  D.sample = rdist(locs.train[,1:2])
  t.sample = rdist(locs.train[,3])
  theta.hat <- c(.9*var(Y.train),mean(D.sample)/4,mean(t.sample)/4,0.1*var(Y.train)) # var,s-range,t-range,nugget
  full.result<-  optim(par=log(theta.hat),fn=negloglik_full_ST,
                       locs=locs.train,y=Y.train,N=p,
                       control=list(trace=FALSE,maxit = 400))
  mdl.geo.fit.avg[i,] <- full.result$par
}

covparam <- exp(apply(mdl.geo.fit.avg,2,mean))

return(covparam)
}


################################################################################
## replicate rows, helps with vector distance
################################################################################
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


################################################################################
## Transform Spatiotemporal data to i.i.d.
################################################################################
transform_iid=function(data,vecchia.approx,covparms,nuggets){
  
  # compute required matrices
  U.obj=createU(vecchia.approx,covparms,nuggets)
  V.ord=U2V(U.obj)
  U.z=U.obj$U[!U.obj$latent,]
  U.y=U.obj$U[U.obj$latent,]
  
  # compute transformed data in parts
  part1.ord=Matrix::crossprod(U.z,data[U.obj$ord.z,])
  temp1=U.y%*%part1.ord
  revord=nrow(temp1):1
  temp2=solve(V.ord,temp1[revord,])
  part2.rev=solve(Matrix::t(V.ord),temp2)
  part2.ord=crossprod(U.y,part2.rev[revord,])
  transform.ord=part1.ord-part2.ord
  
  # return to original ordering
  orig.order=order(U.obj$ord)
  transformed.data=transform.ord[orig.order,]
  return(transform.ord)
  
}

######  GPvecchia local function 
######compute V for posterior inference - needed for transform.iid   #######
U2V=function(U.obj){
  
  U.y=U.obj$U[U.obj$latent,]
  
  if(U.obj$cond.yz=='zy') {
    
    V.ord=revMat(U.y[,U.obj$latent,drop=FALSE])
    
  } else if(U.obj$ord.pred!='obspred'){
    
    W=Matrix::tcrossprod(U.y)
    W.rev=revMat(W)
    V.ord=Matrix::t(Matrix::chol(W.rev))
  } else {  # for obspred ordering
    
    last.obs=max(which(!U.obj$latent))
    latents.before=sum(U.obj$latent[1:last.obs])
    latents.after=sum(U.obj$latent[-(1:last.obs)])
    
    # pred columns are unchanged
    V.pr=revMat(U.y[,(last.obs+1):ncol(U.y),drop=FALSE])
    
    # have to compute cholesky for obs block
    U.oo=U.y[1:latents.before,1:last.obs]
    A=Matrix::tcrossprod(U.oo)
    A.rev=revMat(A)
    V.oor=Matrix::t(Matrix::chol(A.rev))
    
    # combine the blocks into one matrix
    zeromat.sparse=Matrix::sparseMatrix(c(),c(),dims=c(latents.after,latents.before))
    V.or=rbind(zeromat.sparse,V.oor)
    
    V.ord=methods::as(cbind(V.pr,V.or),'dtCMatrix')
    
  }
  
  return(V.ord)
}

###########################################################
## Reverse order of matrix rows,cols
revMat=function(mat) mat[nrow(mat):1,ncol(mat):1,drop=FALSE]





