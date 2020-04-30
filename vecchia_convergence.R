library(MASS)
#library(ggplot2)
#setwd("C:\\Users\\kmessier\\Dropbox\\Papers\\Stats_Paper\\varSelection_kyle")


###### Libraries needed in Vecchia ###########
library(GpGp)
library(Matrix)
library(RcppParallel)
library(parallel)# in order to use detectCores
library(sparseinv)
library(fields)
#
# library(fields) # only needed for pred.cond='independent'
# # Library for SCAD (elastic net)
library(glmnet)
library(ncvreg)
library(geoR)
# install.packages("dequer", repos = "https://cloud.r-project.org")
library(dequer)
# # Libraries for Kriging
#library(geoR)
library(readxl)
library(reshape2)
library(mvtnorm)
# library(optimParallel)
library(scoringRules)

# for (nm in list.files('/home/emt/messierk/R/GPvecchia-0.1.0/R',pattern = "\\.[RrSsQq]$")) {
#   cat(nm,":"); source(file.path('/home/emt/messierk/R/GPvecchia-0.1.0/R',nm)); cat("\n")
# }
# 
# Rcpp::sourceCpp('/home/emt/messierk/R/GPvecchia-0.1.0/src/U_NZentries.cpp')
# Rcpp::sourceCpp('/home/emt/messierk/R/GPvecchia-0.1.0/src/MaxMin.cpp')
# 


for (nm in list.files('C:\\Users\\kmessier\\Dropbox\\Toolboxes\\GPvecchia-master\\R',pattern = "\\.[RrSsQq]$")) {
  cat(nm,":"); source(file.path('C:\\Users\\kmessier\\Dropbox\\Toolboxes\\GPvecchia-master\\R',nm)); cat("\n")
}

Rcpp::sourceCpp('C:\\Users\\kmessier\\Dropbox\\Toolboxes\\GPvecchia-master\\src\\U_NZentries.cpp')
Rcpp::sourceCpp('C:\\Users\\kmessier\\Dropbox\\Toolboxes\\GPvecchia-master\\src\\MaxMin.cpp')
# Rcpp::sourceCpp('C:\\Users\\kmessier\\Dropbox\\Toolboxes\\GPvecchia-master\\src\\fastTree.cpp')
# Rcpp::sourceCpp('C:\\Users\\kmessier\\Dropbox\\Toolboxes\\GPvecchia-master\\src\\ICC.cpp')

## replicate rows, helps with vector distance
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

##### Vecchia functions ################
SCAD.Penalty.Loglike <- function(beta.in,lambda){
  
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


# locs.scaled = cbind(locs[,1]/alpha1,locs[,2]/alpha1,times/alpha2) 
# vecchia.approx = vecchia_specify(locs.scaled, m)

negloglik.vecchia.ST=function(logparms,locs,res,vecchia.approx){
  parms = exp(logparms)
  locs.scaled = cbind(locs[,1]/parms[2], locs[,2]/parms[2], locs[,3]/parms[3])
  vecchia.approx$locsord = locs.scaled[vecchia.approx$ord,]
  -vecchia_likelihood(res, vecchia.approx, c(parms[1],1, 0.5), parms[4] ) 
}  



#### full Kriging loglikelihood [exponential + nugget covariance model]
negloglik.full.ST=function(logparms,locs,y,N){
  parms = exp(logparms)
  locs.scaled = cbind(locs[,1]/parms[2], locs[,2]/parms[2], locs[,3]/parms[3])
  d <- fields::rdist(locs.scaled)
  cov.mat=parms[1]*fields::Exponential(d,range=1)+
    parms[4]*diag(N)
  -mvtnorm::dmvnorm(y,rep(0,N),cov.mat,log=TRUE)
}



### function to transform data to iid using general vecchia
transform.iid=function(data,vecchia.approx,covparms,nuggets){
  
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

#Create a function for Simple Kriging Prediction with a Local neigborhood
Kr.pred  =  function(new_coords,obs_coords,Y_obs,cov.pars,NN){
  
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




##### true parameters
n=500
p=5
beta0=rep(0,p)
theta0=c(9,.3,.5,1) # variance,range,smoothness,nugget


##### simulate data
set.seed(99999)
locs=cbind(runif(n),runif(n)) # random sample on unit square
Sigma.X=exp(-rdist(sample(1:p))/5)
X=mvrnorm(n,rep(0,p),Sigma.X) # correlated predictors (in practice, should include intercept)
Sigma0=theta0[1]*Matern(rdist(locs),range=theta0[2],smoothness=theta0[3])+theta0[4]*diag(n)
epsilon=mvrnorm(1,mu=rep(0,n),Sigma0)
y=X%*%beta0+epsilon


##### convergence for increasing m

# exact transformation
prec.chol=solve(t(chol(Sigma0)))
y.tilde0=prec.chol%*%y
X.tilde0=prec.chol%*%X
beta.hat.exact=solve(crossprod(X.tilde0),crossprod(X.tilde0,y.tilde0))

## vecchia solutions
ms=seq(0,10) #,15,20) #,30,40)
beta.hat.vecchia=matrix(nr=p,nc=length(ms))
for(i.m in 1:length(ms)){
  vecchia.approx=vecchia_specify(locs,ms[i.m])
  transformed.data=transform.iid(cbind(y,X),vecchia.approx,theta0[1:3],theta0[4])
  y.tilde=transformed.data[,1]
  X.tilde=transformed.data[,-1]
  beta.hat.vecchia[,i.m]=as.numeric(solve(crossprod(X.tilde),crossprod(X.tilde,y.tilde)))
}

##### plot of beta estimates as a function of m
pdf(file='plots/vecchiaconv.pdf',width=6.5,height=3.8)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.1,.7)) # bltr
matplot(ms,t(beta.hat.vecchia),type='l',lty=2,lwd=1.5,
        xlab='m',ylab='beta.hat')
abline(h=beta.hat.exact,col=1:p)
legend('bottomright',c('exact','Vecchia'),lty=1:2)
dev.off()

# # plot spatial error field
# quilt.plot(locs[,1],locs[,2],epsilon)
# quilt.plot(locs[,1],locs[,2],X[,1])
# quilt.plot(locs[,1],locs[,2],y)

df <- melt(beta.hat.vecchia)
df.2<- cbind(df,"exact" = rep(beta.hat.exact,11))
df.final <- melt(df.2, id.vars = c("Var1","Var2"))

(p1 <- ggplot(data = df.final,aes(x = Var2,y = value,
                linetype = factor(variable),color = factor(Var1)))+
    geom_line(size=1.25)+labs(y = expression(hat(beta)), x = "m")+
    scale_color_viridis_d(option = "D",aesthetics = c("color"))+
    theme(aspect.ratio = 0.75,legend.position = "bottom")+
    scale_x_continuous(breaks = scales::pretty_breaks(5))+
    scale_linetype_manual(values=c("dashed", "solid"))
  )
