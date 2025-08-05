
# Time cost
TimeCost <- function(time1) {
  timecost <-  proc.time() - time1
  return(timecost)
}

KerNorm<-function(u, h){   
  exp(-u^2/(2*h^2))/h
}
gradKerNorm<-function(u, h){   
  -u*exp(-u^2/(2*h^2))/h^3
}
require(compiler)
KerNormc<- cmpfun(KerNorm)
gradKerNormc<- cmpfun(gradKerNorm)

# pihat estimation 1: using GLM with logit link
EstPI <- function(mydata){
  piestf1 <- glm(delta~mydata$samp1$X, family = binomial, data = mydata$samp1)
  piestf2 <- glm(delta~mydata$samp2$X, family = binomial, data = mydata$samp2)
  pi1hatgamma1 <- piestf1$coefficients
  pi2hatgamma2 <- piestf2$coefficients
  pi1hat <- drop(plogis(piestf1$linear.predictors))## pi1hat
  pi2hat <- drop(plogis(piestf2$linear.predictors))## pi2hat
  
  return(list(pi1hat=pi1hat,pi2hat=pi2hat))
}


# pihat estimation 3: using Random Forest
EstPIrf<- function(mydata){
  y1 <- mydata$samp1$Y
  x1 <- mydata$samp1$X
  delta1 <- mydata$samp1$delta
  
  y2 <- mydata$samp2$Y
  x2 <- mydata$samp2$X
  delta2 <- mydata$samp2$delta
  
  n1 <- length(y1)
  n2 <- length(y2)
  d <- ncol(x1)
  ntree.max <-1000
  tree_vari <- ceiling(d / 2)
  
  if(sum(delta1) == n1 | sum(delta2) == n2 |sum(delta1) == 0 | sum(delta2) == 0){
    pi1hat <- rep(NA, n1)
    pi2hat <- rep(NA, n2)
  }else{
    
    require(randomForest)
    set.seed(666)
    
    piforest1 <- randomForest(x1, y = as.factor(delta1), mtry=tree_vari, ntree=1e+4) #????ɭ??
    
    piforest2 <- randomForest(x2, y = as.factor(delta2), mtry=tree_vari, ntree=1e+4) #????ɭ??
    
    
    pi1hat <- predict(piforest1,x1,"prob")[,2]
    pi2hat <- predict(piforest2,x2,"prob")[,2]
    
    # Extreme value adjustment
    pi1hat[pi1hat-0.001 <= 0] <- 0.001
    pi1hat[pi1hat == 1] <- 0.999
    pi2hat[pi2hat-0.001 <= 0] <- 0.001
    pi2hat[pi2hat == 1] <- 0.999
  }
  return(list(pi1hat=pi1hat,pi2hat=pi2hat))
}

# calculate ispline function of x
cal_bs <- function(x, mm=3, df=6){
  # x is a vector
  require(splines)
  a <- max(abs(x))
  CC <- matrix(0,df-1,df)   
  Boundary.knots <- range(-a, a)
  order <- mm + 1
  nIknots <- df-order
  knots <- seq.int(from =(-a), to=a, length.out = nIknots + 2L)[-c(1L, nIknots + 2L)]  
  BS1 <- bs(x,degree=mm,knots=knots,Boundary.knots=Boundary.knots, intercept =T)  
  return(list(BS=BS1))
}

# simplify the xls
ClearXls <- function(x,mydata){
  y1 <- mydata$samp1$Y
  y2 <- mydata$samp2$Y
  n1 <- length(y1)
  n2 <- length(y2)
  
  dy <- outer(c(y1),c(y2),"-")
  dy.max <- max(dy)
  dy.min <- min(dy)
  
  xnew <- x[-1]# drop x=0
  xn0 <- length(which(xnew<dy.min))
  xn1 <- length(which(xnew>dy.max))
  xbetween <- xnew[which(dy.min<=xnew& xnew<=dy.max)]
  
  return(list(xn0=xn0,xcompu=c(0,xbetween),xn1=xn1))
}

# Fx estimation
EstFx <- function(xlist, mydata, thetahat, hatpi, df=6, flag=1, hm){
  ## thetahat is a matrix
  X1 <- mydata$samp1$X
  X2 <- mydata$samp2$X
  Y1 <- mydata$samp1$Y
  Y2 <- mydata$samp2$Y
  delta1 <- mydata$samp1$delta
  delta2 <- mydata$samp2$delta
  
  p1=p2=ncol(X1)
  n <- length(Y1)
  m <- length(Y2)
  
  pi1 <- hatpi$pi1hat
  pi2 <- hatpi$pi2hat
  
  theta1=thetahat[1:(p1*df)]
  theta2=thetahat[(p1*df+1):(p1*df+p2*df)]
  if(flag==1){
    indtheta=seq(df,p1*df,by=df)
    beta1=theta1[indtheta]
    beta2=theta2[indtheta]
    sigma_betahatx=sqrt(t(beta1)%*%cov(X1)%*%beta1 + t(beta2)%*%cov(X2)%*%beta2)
    hm=c(sigma_betahatx)*(n+m)^{-1/3}  
  }
  
  # data cleaning
  delta12.1=c(outer(c(delta1), c(delta2)))
  Y12.1=c(outer(c(Y1),c(Y2),"-"))
  Y12=Y12.1[delta12.1==1]
  pi12=outer(c(pi1), c(pi2))[delta12.1==1]
  
  bs.beta <- cal_bs(xlist)$BS
  len.n=length(xlist)
  betai1=matrix(0,len.n,p1)
  for(i1 in 1:p1){
    if(i1==1) starti=1 else starti=(i1-1)*df+1
    endi=i1*df
    thetai1=theta1[starti:endi]
    betai1[,i1]=crossprod(t(bs.beta),thetai1)
  }
  betai2=matrix(0,len.n,p2)
  for(i2 in 1:p2){
    if(i2==1) starti=1 else starti=(i2-1)*df+1
    endi=i2*df
    thetai2=theta2[starti:endi]
    betai2[,i2]=crossprod(t(bs.beta),thetai2)
  }
  Sx1.1=tcrossprod(X1,betai1)   ### X_1'\beta_1(x); dim=n*len.n
  Sx2.1=tcrossprod(X2,betai2)   ### X_2'\beta_2(x); dim=m*len.n
  Sx.1=matrix(0,m*n,len.n)
  for(i in 1:n){
    if(i==1) starti=1 else starti=(i-1)*m+1
    endi=i*m
    Sx.1[starti:endi,]=t(t(Sx1.1)+Sx2.1[i,])
  }
  Sx=Sx.1[delta12.1==1,]   ### DIM=(Sx1i+SXj2)*point(x)
  
  XX <- rbind(X1,X2)
  SSx1.1=tcrossprod(XX,betai1)   
  SSx2.1=tcrossprod(XX,betai2)  
  SSx=matrix(0,(m+n)*(m+n),len.n)
  for(i in 1:(m+n)){
    if(i==1) starti=1 else starti=(i-1)*(m+n)+1
    endi=i*(m+n)
    SSx[starti:endi,]=t(t(SSx1.1)+SSx2.1[i,])
  }
  
  Fn=rep(0,len.n)
  for(ix in 1:len.n){
    x=xlist[ix]
    indY <- ifelse(Y12<=x, 1, 0)
    Sx.temp=Sx[,ix]
    SSx.temp=SSx[,ix]
    matSx.temp <- outer(Sx.temp, SSx.temp, "-")
    xKmat <- KerNormc(matSx.temp,hm)  ## K(Sij(x)-s) for fixed x
    m1val <- colSums(1/pi12*indY*xKmat)
    m2val <- colSums(1/pi12*xKmat)
    m2val[m2val==0]=1
    Fn[ix] <- mean(m1val/m2val)
  }
  return(Fn)
}

EstBetaRS2 <- function(mydata, xlist,hatpi){
  pi1hat <- hatpi$pi1hat
  pi2hat <- hatpi$pi2hat
  
  x1 <- as.matrix(mydata$samp1$X)
  x2 <- as.matrix(mydata$samp2$X)
  y1 <- mydata$samp1$Y
  y2 <- mydata$samp2$Y
  delta1 <- mydata$samp1$delta
  delta2 <- mydata$samp2$delta
  n1 <- length(y1)
  n2 <- length(y2)
  d <- ncol(x1)
  
  
  deltapi <- c(outer(c(delta1), c(delta2)) / outer(c(pi1hat), c(pi2hat)))
  dy <- c(outer(c(y1), c(y2), "-"))
  nonzero.index <- which(deltapi!=0)
  deltapi.vec <- deltapi[nonzero.index]
  dy <- dy[nonzero.index]
  
  YY0 <- ifelse(dy <= quantile(dy, 0.5), 1, 0)
  XX.init <- cbind(kronecker(rep(1, n2), x1), kronecker(x2, rep(1, n1)))## D11,D21,D31,...,D12,D22,D32,...
  XX <- XX.init[nonzero.index,]
  
  betainit1 <- glm(YY0~XX+0, family =binomial, weights=1/deltapi.vec)
  if(betainit1$converged){
    betainit=betainit1$coefficients
    betainit[is.na(betainit)]=0
    betainit[is.na(betainit)]=0
  }  else betainit=rep(NA,2*d)
  
  thetahat=numeric(0)
  if(!anyNA(betainit)){
    
    betahat.mat <- matrix(NA, 2 * d, length(xlist))
    
    for (i in 1:length(xlist)) {
      YY <- ifelse(dy <= xlist[i], 1, 0)
      dataloss <- list(y=YY, x=XX)
      
      require(splines)
      lossFx <- function(beta,dataloss,dspline=6){
        y <- dataloss$y
        x <- dataloss$x
        s <- bs(x%*%beta, df = dspline)
        yhat <- fitted(lm(y~s-1,weights=1/deltapi.vec),s)
        lval <- sum((y-yhat)^2/deltapi.vec,na.rm = T)
        return(lval)
      }
      
      require(nlme)
      betahat.mat[,i] <- nlm(lossFx, betainit, data=dataloss)$estimate
    }
    sx <- bs(xlist, df =6)
    for(j in 1:(2*d)){
      thetahat=c(thetahat,lm(betahat.mat[j,]~sx+0)$coef)
    }
    thetahat[is.na(thetahat)]=0
    thetahat[is.nan(thetahat)]=0
  }else{
    thetahat <- rep(1,2*d*df)
  }
  return(thetahat)
}

