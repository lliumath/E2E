
# space cleaning
ls()
remove(list = ls())

# files sourcing
# ------------------------------------
source("data_generation.R")
source("combined_methods.R")
source("proposed_betax.R")


#----------gdata1-6-------------------------------------
nlist <- c(50,100)
plist <- c(5)
ssize <- 400
kerf <- "norm"
misstype="MAR1"

xgrid.mun <- 50
xlist=c(0, seq(-22,22,length.out=xgrid.mun)) 
gcase <- "gdata4"
yroom <- "room4"
ws.init <- paste("/home/", yroom, sep = "")
mlist= c("fxbetaxpro","fxbetaxprorf")
mi.len=length(mlist)


for (mi in 1:mi.len) {
  dir.create(file.path(ws.init,mlist[mi]))

  tclist <- NA
  tl <- 1
  for (n in nlist) {
    dir.create(file.path(paste(ws.init,"/",mlist[mi],sep = ""),n))
    setwd(paste(ws.init,"/",mlist[mi],"/",n,sep = ""))
    for (p in plist) {
       vls=list(n1=n, n2=n, p=p, misstype=misstype, gcase=gcase, xls=xlist)
       cores = 30

  time1 <- proc.time()
  
  require(doSNOW)
  require(tcltk)
  require(doParallel)
  cls <- makeSOCKcluster(cores, revtunnel = TRUE,
                         outfile = paste("/home/liuli/ll/", yroom, "/out.log", sep = ""), verbose = TRUE)
  
  registerDoSNOW(cls)
  
  pb <- txtProgressBar(min=1, max=ssize, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  cat(paste("\033[31m Serial Number:\033[0m",tl),
      paste("\033[32m method:\033[0m",mlist[mi]),
      paste("\033[33m n:\033[0m",n),
      paste("\033[36m d:\033[0m",p),
      "\n")
  ffnx1 <- foreach(i=1:ssize, .options.snow=opts,.combine = "cbind")%dopar% methodTry(mlist[mi],i,vls)

  if(tl%/%10==0){
    fname <- paste("0",tl,"_",ssize,"_",mlist[mi],"_",kerf,"_",n,"_",p,"_mar1.txt",sep="")
  }else{
    fname <- paste(tl,"_",ssize,"_",mlist[mi],"_",kerf,"_",n,"_",p,"_mar1.txt",sep="")
  }
  
  write.table(ffnx1,paste(ws.init,"/",mlist[mi],"/",fname,sep = ""))
  
  close(pb)
  stopCluster(cls)
  
  rm("ffnx1","cls")
  
  tclist[tl] <- TimeCost(time1)[3]
  tl <- tl+1
   }
  }

  tcinfor <- list(tclist=tclist,cores=cores,ssize=ssize)
  write.csv(tcinfor,paste(mlist[mi],"_time_cost.csv",sep=""))
}
