
#///////////////////////////////////////////////////////#
AR <- function(n,p,rho){
  wn.sd <- sqrt(1 - rho^2)
  wn.mat <- matrix(rnorm(n*p), nrow=n, ncol=p) * wn.sd  ## wn means white noise
  Xmat <- matrix(0, nrow=n, ncol=p)
  Xmat[,1] <- rnorm(n)
  for(j in 2:p)
    Xmat[,j] <- rho * Xmat[,j-1] + wn.mat[,j]
  Xmat <- matrix(pmax(-1,pmin(1,Xmat)),n,p)
  return(Xmat)
}

#---------------------------------------------------------
#* gdata1:
#* y1: Normal distribution with mean abs(x1*beta1)
#* y2: Exponential distribution with mean abs(x2*beta2)
#* MISSING: LOGIT
#* rho1=0.2; rho2=0.8
#---------------------------------------------------------
gdata1 <- function(i, n,m, p,misstype="MAR1"){
  #* settings: 
  #* rho1, rho2, 
  #* beta1.true, beta2.true, 
  #* gamma1, gamma2
  set.seed(6000+80*i)
  
  rho1 <- 0.4
  beta1 <- rep(1,p)
  beta2 <- 0.3*rep(1,p)
  X1.1=X2.1=X=AR(n,p,rho1)
  
  X.2 <- cbind(rep(1,n), X)
  gamma=c(-0.2,rep(1/sqrt(p),p))
  missp.1 <- apply(t(X.2)*gamma,2,sum)
  missp <- 1-1/(1+exp(missp.1))
  delta <- rbinom(n,1,missp)
  
  X1 <- c(crossprod(t(X1.1),beta1))
  X2 <- c(crossprod(t(X2.1),beta2))
  
  Y1 <- rnorm(n,X1,1) ## Y1 with sample size n1
  Y2 <- rnorm(n,X2,1) ## Y1 with sample size n1
  
  samp1 <- list(X=X1.1,Y=Y1,delta=delta)
  samp2 <- list(X=X2.1,Y=Y2,delta=1-delta)
  misrate1 <- 1-mean(delta)
  misrate2 <- 1-misrate1
  return(list(samp1=samp1, samp2=samp2, mrate1=misrate1, mrate2=misrate2))
}

gdata2 <- function(i, n,m, p,misstype="MAR1"){
  #* settings: 
  #* rho1, rho2, 
  #* beta1.true, beta2.true, 
  #* gamma1, gamma2
  set.seed(6000+80*i)
  
  rho1 <- 0.4
  beta1 <- rep(1,p)
  beta2 <- 0.3*rep(1,p)
  X1.1=X2.1=X=AR(n,p,rho1)
  
  X.2 <- cbind(rep(1,n), X)
  gamma=c(-0.2,rep(1/sqrt(p),p))
  missp.1 <- apply(t(X.2)*gamma,2,sum)
  missp <- 1-1/(1+exp(missp.1))
  delta <- rbinom(n,1,missp)
  
  X1 <- c(crossprod(t(X),beta1))
  X2 <- c(crossprod(t(abs(X)),beta2))
  X2[X2==0]=mean(X2)
  
  Y1 <- abs(rnorm(n,X1,1)) ## Y1 with sample size n1
  Y2 <- rexp(m,X2*0.3)        ## Y2 with sample size n2
  
  samp1 <- list(X=X1.1,Y=Y1,delta=delta)
  samp2 <- list(X=X2.1,Y=Y2,delta=1-delta)
  misrate1 <- 1-mean(delta)
  misrate2 <- 1-misrate1
  return(list(samp1=samp1, samp2=samp2, mrate1=misrate1, mrate2=misrate2))
}

gdata3 <- function(i, n,m, p,misstype="MAR1"){
  #* settings: 
  #* rho1, rho2, 
  #* beta1.true, beta2.true, 
  #* gamma1, gamma2
  set.seed(6000+80*i)
  
  rho1 <- 0.4
  beta1 <- rep(1,p)
  beta2 <- 0.3*rep(1,p)
  X1.1=X2.1=X=AR(n,p,rho1)
  
  X.2 <- cbind(rep(1,n), X, X[,1]^2, X[,p]^2)
  gamma=c(-0.1,rep(1/sqrt(p),p+2))
  missp.1 <- apply(t(X.2)*gamma,2,sum)
  missp <- 1-1/(1+exp(missp.1))
  delta <- rbinom(n,1,missp)
  
  X1 <- c(crossprod(t(X1.1),beta1))
  X2 <- c(crossprod(t(X2.1),beta2))
  
  Y1 <- rnorm(n,X1,1) ## Y1 with sample size n1
  Y2 <- rnorm(n,X2,1) ## Y1 with sample size n1
  
  samp1 <- list(X=X1.1,Y=Y1,delta=delta)
  samp2 <- list(X=X2.1,Y=Y2,delta=1-delta)
  misrate1 <- 1-mean(delta)
  misrate2 <- 1-misrate1
  return(list(samp1=samp1, samp2=samp2, mrate1=misrate1, mrate2=misrate2))
}

gdata4 <- function(i, n, m, p,misstype="MAR1"){
  #* settings: 
  #* rho1, rho2, 
  #* beta1.true, beta2.true, 
  #* gamma1, gamma2
  set.seed(6000+80*i)
  
  rho1 <- 0.4
  beta1 <- rep(1,p)
  beta2 <- 0.3*rep(1,p)
  X1.1=X2.1=X=AR(n,p,rho1)
  
  X.2 <- cbind(rep(1,n), X, X[,1]^2, X[,p]^2)
  gamma=c(-0.1,rep(1/sqrt(p),p+2))
  missp.1 <- apply(t(X.2)*gamma,2,sum)
  missp <- 1-1/(1+exp(missp.1))
  delta <- rbinom(n,1,missp)
  
  X1 <- c(crossprod(t(X),beta1))
  X2 <- c(crossprod(t(abs(X)),beta2))
  X2[X2==0]=mean(X2)
  
  Y1 <- abs(rnorm(n,X1,1)) ## Y1 with sample size n1
  Y2 <- rexp(m,X2*0.3)        ## Y2 with sample size n2
  
  samp1 <- list(X=X1.1,Y=Y1,delta=delta)
  samp2 <- list(X=X2.1,Y=Y2,delta=1-delta)
  misrate1 <- 1-mean(delta)
  misrate2 <- 1-misrate1
  return(list(samp1=samp1, samp2=samp2, mrate1=misrate1, mrate2=misrate2))
}

gdata5 <- function(i, n, m, p,misstype="MAR1"){
  #* settings: 
  #* rho1, rho2, 
  #* beta1.true, beta2.true, 
  #* gamma1, gamma2
  set.seed(6000+80*i)
  
  rho1 <- 0.4
  rho=1
  beta1 <- rep(1,p)
  beta2 <- 0.3*rep(1,p)
  Z=AR(n,p,rho1)
  n1=floor(n/2)
  X1.1=Z[1:n1,]-rho #*exp(Z)
  X2.1=Z[(n1+1):n,]+rho #*exp(Z)
  X=rbind(X1.1,X2.1)
  
  X.2 <- cbind(rep(1,n), X)
  #  gamma=c(1.6,rep(1/sqrt(p),p+2))
  gamma=c(-1,rep(1/sqrt(p),p))
  missp.1 <- apply(t(X.2)*gamma,2,sum)
  missp <- 1-1/(1+exp(missp.1))
  delta <- rbinom(n,1,missp)
  
  X1 <- c(crossprod(t(X),beta1))
  X2 <- c(crossprod(t(X),beta2))
  
  Y1 <- rnorm(n,X1,1) ## Y1 with sample size n1
  Y2 <- rnorm(n,X2,1) ## Y1 with sample size n1
  
  samp1 <- list(X=X,Y=Y1,delta=delta)
  samp2 <- list(X=X,Y=Y2,delta=1-delta)
  misrate1 <- 1-mean(delta)
  misrate2 <- 1-misrate1
  return(list(samp1=samp1, samp2=samp2, mrate1=misrate1, mrate2=misrate2))
}

gdata6 <- function(i, n, m, p,misstype="MAR1"){
  #* settings: 
  #* rho1, rho2, 
  #* beta1.true, beta2.true, 
  #* gamma1, gamma2
  set.seed(6000+80*i)
  
  rho1 <- 0.4
  rho=1
  beta1 <- rep(1,p)
  beta2 <- 0.3*rep(1,p)
  Z=AR(n,p,rho1)
  n1=floor(n/2)
  X1.1=Z[1:n1,]-rho #*exp(Z)
  X2.1=Z[(n1+1):n,]+rho #*exp(Z)
  X=rbind(X1.1,X2.1)
  
  X.2 <- cbind(rep(1,n), X)
  #  gamma=c(1.6,rep(1/sqrt(p),p+2))
  gamma=c(-1,rep(1/sqrt(p),p))
  missp.1 <- apply(t(X.2)*gamma,2,sum)
  missp <- 1-1/(1+exp(missp.1))
  delta <- rbinom(n,1,missp)
  
  X1 <- c(crossprod(t(X),beta1))
  X2 <- c(crossprod(t(abs(X)),beta2))
  X2[X2==0]=mean(X2)
  
  Y1 <- abs(rnorm(n,X1,1)) ## Y1 with sample size n1
  Y2 <- rexp(m,X2*0.3)        ## Y2 with sample size n2
  
  samp1 <- list(X=X,Y=Y1,delta=delta)
  samp2 <- list(X=X,Y=Y2,delta=1-delta)
  misrate1 <- 1-mean(delta)
  misrate2 <- 1-misrate1
  return(list(samp1=samp1, samp2=samp2, mrate1=misrate1, mrate2=misrate2))
}

gdata7<- function(i, n,m, p,misstype="MAR1"){
  #* settings: 
  #* rho1, rho2, 
  #* beta1.true, beta2.true, 
  #* gamma1, gamma2
  set.seed(6000+80*i)
  
  rho1 <- 0.4
  rho=1
  beta1 <- rep(1,p)
  beta2 <- 0.3*rep(1,p)
  Z=AR(n,p,rho1)
  n1=floor(n/2)
  X1.1=Z[1:n1,]-rho #*exp(Z)
  X2.1=Z[(n1+1):n,]+rho #*exp(Z)
  X=rbind(X1.1,X2.1)
  
  X.2 <- cbind(rep(1,n), X, X[,1]^2, X[,p]^2)
#  gamma=c(1.6,rep(1/sqrt(p),p+2))
  gamma=c(-1,rep(1/sqrt(p),p+2))
  missp.1 <- apply(t(X.2)*gamma,2,sum)
  missp <- 1-1/(1+exp(missp.1))
  delta <- rbinom(n,1,missp)
  
  X1 <- c(crossprod(t(X),beta1))
  X2 <- c(crossprod(t(X),beta2))
  
  Y1 <- rnorm(n,X1,1) ## Y1 with sample size n1
  Y2 <- rnorm(n,X2,1) ## Y1 with sample size n1
  
  samp1 <- list(X=X,Y=Y1,delta=delta)
  samp2 <- list(X=X,Y=Y2,delta=1-delta)
  misrate1 <- 1-mean(delta)
  misrate2 <- 1-misrate1
  return(list(samp1=samp1, samp2=samp2, mrate1=misrate1, mrate2=misrate2))
}

gdata8<- function(i, n,m, p,misstype="MAR1"){
  #* settings: 
  #* rho1, rho2, 
  #* beta1.true, beta2.true, 
  #* gamma1, gamma2
  set.seed(6000+80*i)
  
  rho1 <- 0.4
  rho=1
  beta1 <- rep(1,p)
  beta2 <- 0.3*rep(1,p)
  Z=AR(n,p,rho1)
  n1=floor(n/2)
  X1.1=Z[1:n1,]-rho #*exp(Z)
  X2.1=Z[(n1+1):n,]+rho #*exp(Z)
  X=rbind(X1.1,X2.1)
  
  X.2 <- cbind(rep(1,n), X, X[,1]^2, X[,p]^2)
  #  gamma=c(1.6,rep(1/sqrt(p),p+2))
  gamma=c(-1,rep(1/sqrt(p),p+2))
  missp.1 <- apply(t(X.2)*gamma,2,sum)
  missp <- 1-1/(1+exp(missp.1))
  delta <- rbinom(n,1,missp)
  
  X1 <- c(crossprod(t(X),beta1))
  X2 <- c(crossprod(t(abs(X)),beta2))
  X2[X2==0]=mean(X2)
  
  Y1 <- abs(rnorm(n,X1,1)) ## Y1 with sample size n1
  Y2 <- rexp(m,X2*0.3)        ## Y2 with sample size n2
  
  samp1 <- list(X=X,Y=Y1,delta=delta)
  samp2 <- list(X=X,Y=Y2,delta=1-delta)
  misrate1 <- 1-mean(delta)
  misrate2 <- 1-misrate1
  return(list(samp1=samp1, samp2=samp2, mrate1=misrate1, mrate2=misrate2))
}

gdata7<- function(i, n,m, p,misstype="MAR1"){
  #* settings: 
  #* rho1, rho2, 
  #* beta1.true, beta2.true, 
  #* gamma1, gamma2
  set.seed(6000+80*i)
  
  rho1 <- 0.4
  rho=1
  beta1 <- rep(1,p)
  beta2 <- 0.3*rep(1,p)
  Z=AR(n,p,rho1)
  n1=floor(n/2)
  X1.1=Z[1:n1,]-rho #*exp(Z)
  X2.1=Z[(n1+1):n,]+rho #*exp(Z)
  X=rbind(X1.1,X2.1)
  
  X.2.1=exp(X[,1]/4)  
  X.2.2=X[,2]/(1+exp(2*X[,1]))
  X.2.3=(X[,2]*X[,3]-1)^3
  X.2.4=(X[,3]+X[,4]+X[,5])^2/5
  X.2 <- cbind(rep(1,n), X.2.1, X.2.2, X.2.3, X.2.4)
  #  gamma=c(1.6,rep(1/sqrt(p),p+2))
  gamma=c(-1.5,rep(1/sqrt(p),p-1))
  missp.1 <- apply(t(X.2)*gamma,2,sum)
  missp <- 1-1/(1+exp(missp.1))
  delta <- rbinom(n,1,missp)
  
  X1 <- c(crossprod(t(X),beta1))
  X2 <- c(crossprod(t(X),beta2))
  
  Y1 <- rnorm(n,X1,1) ## Y1 with sample size n1
  Y2 <- rnorm(n,X2,1) ## Y1 with sample size n1
  
  samp1 <- list(X=X,Y=Y1,delta=delta)
  samp2 <- list(X=X,Y=Y2,delta=1-delta)
  misrate1 <- 1-mean(delta)
  misrate2 <- 1-misrate1
  return(list(samp1=samp1, samp2=samp2, mrate1=misrate1, mrate2=misrate2))
}


gdata8<- function(i, n,m, p,misstype="MAR1"){
  #* settings: 
  #* rho1, rho2, 
  #* beta1.true, beta2.true, 
  #* gamma1, gamma2
  set.seed(6000+80*i)
  
  rho1 <- 0.4
  rho=1
  beta1 <- rep(1,p)
  beta2 <- 0.3*rep(1,p)
  Z=AR(n,p,rho1)
  n1=floor(n/2)
  X1.1=Z[1:n1,]-rho #*exp(Z)
  X2.1=Z[(n1+1):n,]+rho #*exp(Z)
  X=rbind(X1.1,X2.1)
  
  X.2.1=exp(X[,1]/4)  
  X.2.2=X[,2]/(1+exp(2*X[,1]))
  X.2.3=(X[,2]*X[,3]-1)^3
  X.2.4=(X[,3]+X[,4]+X[,5])^2/5
  X.2 <- cbind(rep(1,n), X.2.1, X.2.2, X.2.3, X.2.4)
  #  gamma=c(1.6,rep(1/sqrt(p),p+2))
  gamma=c(-1.5,rep(1/sqrt(p),p-1))
  missp.1 <- apply(t(X.2)*gamma,2,sum)
  missp <- 1-1/(1+exp(missp.1))
  delta <- rbinom(n,1,missp)
  
  X1 <- c(crossprod(t(X),beta1))
  X2 <- c(crossprod(t(abs(X)),beta2))
  X2[X2==0]=mean(X2)
  
  Y1 <- abs(rnorm(n,X1,1)) ## Y1 with sample size n1
  Y2 <- rexp(m,X2*0.3)        ## Y2 with sample size n2
  
  samp1 <- list(X=X,Y=Y1,delta=delta)
  samp2 <- list(X=X,Y=Y2,delta=1-delta)
  misrate1 <- 1-mean(delta)
  misrate2 <- 1-misrate1
  return(list(samp1=samp1, samp2=samp2, mrate1=misrate1, mrate2=misrate2))
}
