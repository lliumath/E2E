
# Fx Method-"fxbetaxpro"(FxBetax-proposed hm=0.3)
wfun0 <- function(i,vls){
  flag=0
  hm=0.3
  if(vls$gcase=="gdata1"){
    datai <- gdata1(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata2"){
    datai <- gdata2(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata3"){
    datai <- gdata3(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata4"){
    datai <- gdata4(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata5"){
    datai <- gdata5(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata6"){
    datai <- gdata6(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata7"){
    datai <- gdata7(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata8"){
    datai <- gdata8(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  pihat <- EstPI(datai)
  pi1hat <- pihat$pi1hat
  pi2hat <- pihat$pi2hat
  
  if(!is.na(sum(pi1hat))&!is.na(sum(pi2hat))){
    cxls <- ClearXls(vls$xls,datai)
    xbeta.ls <- c(0,quantile(cxls$xcompu[-1],seq(0,1,length.out = 30)))
    thetahat <- EstBetaRS2(datai,xbeta.ls,pihat)
    fbetween <- EstFx(cxls$xcompu, datai, thetahat, pihat, df=6, flag, hm)
    ffn1 <- c(fbetween[1],rep(0,cxls$xn0),fbetween[-1],rep(1,cxls$xn1))
   }

  if(i%/%10==0){
    fname <- paste("000",i,".txt",sep="")
  }else if(i%/%100==0){
    fname <- paste("00",i,".txt",sep="")
  }else if(i%/%1000==0){
    fname <- paste("0",i,".txt",sep="")
  }else{
    fname <- paste(i,".txt",sep="")
  }
  
  write.table(ffn1,fname)
  return(ffn1)
}


# Fx Method-"fxbetaxprorf6"
wfun11 <- function(i,vls){
  flag=1
  hm=0.3
  if(vls$gcase=="gdata1"){
    datai <- gdata1(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata2"){
    datai <- gdata2(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata3"){
    datai <- gdata3(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata4"){
    datai <- gdata4(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata5"){
    datai <- gdata5(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata6"){
    datai <- gdata6(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata7"){
    datai <- gdata7(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  if(vls$gcase=="gdata8"){
    datai <- gdata8(i, vls$n1, vls$n2, vls$p, vls$misstype)
  }
  pihat <- EstPIrf(i,datai)

  pi1hat <- pihat$pi1hat
  pi2hat <- pihat$pi2hat
  
  if(!is.na(sum(pi1hat))&!is.na(sum(pi2hat))){
    cxls <- ClearXls(vls$xls,datai)
    xbeta.ls <- c(0,quantile(cxls$xcompu[-1],seq(0,1,length.out = 30)))
    thetahat <- EstBetaRS2(datai,xbeta.ls,pihat)
    fbetween <- EstFx(cxls$xcompu, datai, thetahat, pihat, df=6, flag, hm)
    ffn1 <- c(fbetween[1],rep(0,cxls$xn0),fbetween[-1],rep(1,cxls$xn1))
  }

  if(i%/%10==0){
    fname <- paste("000",i,".txt",sep="")
  }else if(i%/%100==0){
    fname <- paste("00",i,".txt",sep="")
  }else if(i%/%1000==0){
    fname <- paste("0",i,".txt",sep="")
  }else{
    fname <- paste(i,".txt",sep="")
  }
  
  write.table(ffn1,fname)
   
  return(ffn1)
}


# Methods collection
methodTry <- function(mname,i,vls){
  
  if(mname=="fxbetaxpro"){
    mname <- "fxbetaxpro"
    reval <- wfun0(i,vls)
  }
  
  if(mname=="fxbetaxprorf"){
    mname <- "fxbetaxprorf"
    reval <- wfun11(i,vls)
  }
  return(reval)
}


