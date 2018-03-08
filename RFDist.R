# RF Dissimilarities
library(MASS)
library(cluster)
library(survival)
library(randomForest)
library(Hmisc)

collect.garbage <- function(){
  ## The following function collects garbage until the memory is clean.
  ## Usage: 1. immediately call this function after you call a function or
  ##        2. rm()
  while (gc()[2,4] != gc()[2,4]){}
}



RFdist <- function(datRF, mtry1, no.tree, no.forest, addcl1=T, addcl2=T, imp=T, oob.prox1=T, proxConver=F) {

# dealing with unlabeled data  
# added by randomly sampling from the product of empirical marginal dist of variables  
  synthetic1 <- function(dat) {
    sample1 <- function(X)   { sample(X, replace=T) } 
    g1      <- function(dat) { apply(dat,2,sample1) }
    nrow1 <- dim(dat)[[1]];
    yy <- rep(c(1,2),c(nrow1,nrow1) );
    data.frame(cbind(yy,rbind(dat,data.frame(g1(dat)))))
  }

# added by sampling from the hyper-rectangle that contains the obs data(uniform dist with max and min of corresponding obs data)  
  synthetic2 <- function(dat) {
    sample2 <- function(X)   { runif(length(X), min=min(X), max =max(X)) }
    g2      <- function(dat) { apply(dat,2,sample2) }
    nrow1 <- dim(dat)[[1]];
    yy <- rep(c(1,2),c(nrow1,nrow1) );
    data.frame(cbind(yy,rbind(dat,data.frame(g2(dat)))))
  }
 
  synthetic3 <- function(dat) {
    sample3 <- function(X)   { rnorm(length(X), mean=mean(X), sd = var(X)) }
    g3      <- function(dat) { apply(dat,2,sample3) }
    nrow1 <- dim(dat)[[1]];
    yy <- rep(c(1,2),c(nrow1,nrow1) );
    data.frame(cbind(yy,rbind(dat,data.frame(g3(dat)))))
  }
   
  cleandist <- function(x) { 
    x1 <- as.dist(x)
    x1[x1<=0] <- 0.0000000001
    as.matrix(x1)
  }
  
  nrow1 <- dim(datRF)[[1]]
  ncol1 <- dim(datRF)[[2]]
  RFproxAddcl1 <- matrix(0,nrow=nrow1,ncol=nrow1)
  RFproxAddcl2 <- matrix(0,nrow=nrow1,ncol=nrow1)
  RFprox1Conver <- cbind(1:no.forest,matrix(0,(no.forest),3))
  RFprox2Conver <- cbind(1:no.forest,matrix(0,(no.forest),3))
  RFimportance1 <- matrix(0, nrow=ncol1, ncol=4)
  RFimportance2 <- matrix(0, nrow=ncol1, ncol=4)
  RFerrrate1 <- 0
  RFerrrate2 <- 0
  rep1 <- rep(666,2*nrow1) 
  
  #addcl1
  if (addcl1) {
    for (i in c(0:no.forest)) { 
      index1 <- sample(c(1:(2*nrow1))) 
      rep1[index1] <-  c(1:(2*nrow1)) 
      datRFsyn <- synthetic1(datRF)[index1,] 
      yy <- datRFsyn[,1] 
      RF1 <- randomForest(factor(yy)~.,data=datRFsyn[,-1], ntree=no.tree, oob.prox=oob.prox1, proximity=TRUE,do.trace=F,mtry=mtry1,importance=imp) 
      collect.garbage()
      RF1prox <- RF1$proximity[rep1,rep1]
      if (i > 0) { 
        if (i > 1){
          xx <- ((RFproxAddcl1 + (RF1prox[c(1:nrow1),c(1:nrow1)]))/i) - (RFproxAddcl1/(i-1))
          yy <- mean( c(as.dist((RFproxAddcl1 + (RF1prox[c(1:nrow1),c(1:nrow1)]))/i))) 
          RFprox1Conver[i,2] <- max(abs(c(as.dist(xx))))
          RFprox1Conver[i,3] <- mean((c(as.dist(xx)))^2)
          RFprox1Conver[i,4] <- yy
        }
        RFproxAddcl1 <- RFproxAddcl1 + (RF1prox[c(1:nrow1),c(1:nrow1)]) 
        if(imp) { RFimportance1 <- RFimportance1+ 1/no.forest*(RF1$importance) }
        RFerrrate1 <- RFerrrate1+ 1/no.forest*(RF1$err.rate[no.tree])
      }
    }
  }
  
  # addcl2
  if (addcl2) { 
    for (i in c(0:no.forest)) {
      index1 <- sample(c(1:(2*nrow1))) 
      rep1[index1] <-  c(1:(2*nrow1)) 
      datRFsyn <- synthetic2(datRF)[index1,] 
      yy <- datRFsyn[,1] 
      RF2 <- randomForest(factor(yy)~.,data=datRFsyn[,-1], ntree=no.tree, oob.prox=oob.prox1, proximity=TRUE,do.trace=F,mtry=mtry1,importance=imp) 
      collect.garbage()
      RF2prox <- RF2$proximity[rep1,rep1]
      if (i > 0) { 
        if (i > 1){
          xx <- ((RFproxAddcl2 + (RF2prox[c(1:nrow1),c(1:nrow1)]))/i) - (RFproxAddcl2/(i-1))
          yy <- mean( c(as.dist((RFproxAddcl2 + (RF2prox[c(1:nrow1),c(1:nrow1)]))/i))) 
          RFprox2Conver[i,2] <- max(abs(c(as.dist(xx))))
          RFprox2Conver[i,3] <- mean((c(as.dist(xx)))^2)
          RFprox2Conver[i,4] <- yy
        }
        RFproxAddcl2 <- RFproxAddcl2 + (RF2prox[c(1:nrow1),c(1:nrow1)]) 
        if(imp) { RFimportance2 <- RFimportance2+ 1/no.forest*(RF2$importance)}
        RFerrrate2 <- RFerrrate2+ 1/no.forest*(RF2$err.rate[no.tree])
      }
    }
  }
  
#  cl1:  addcl1 distance (sqrt(1-RF.proxAddcl1))
#  cl2:  addcl2 distance (sqrt(1-RF.proxAddcl2))
  
  distRFAddcl1 <- cleandist(sqrt(1-RFproxAddcl1/no.forest))
  distRFAddcl2 <- cleandist(sqrt(1-RFproxAddcl2/no.forest))
  
  distRF <- list(cl1=NULL, err1=NULL, imp1=NULL, prox1Conver=NULL, 
                 cl2=NULL, err2=NULL, imp2=NULL, prox2Conver=NULL)
  if(addcl1) {
    distRF$cl1 <- distRFAddcl1
    distRF$err1 <- RFerrrate1
    if(imp) distRF$imp1 <- RFimportance1 
    if(proxConver) distRF$prox1Conver <- RFprox1Conver
  }
  if(addcl2) {
    distRF$cl2 <- distRFAddcl2
    distRF$err2 <- RFerrrate2
    if(imp) distRF$imp2 <- RFimportance2
    if(proxConver) distRF$prox2Conver <- RFprox2Conver
  } 
  distRF
}

