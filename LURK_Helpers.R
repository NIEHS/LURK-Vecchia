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
    print(i)
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
