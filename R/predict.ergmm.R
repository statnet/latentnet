predict.ergmm<-function(ergmm.fit,theta="mkl"){
  if(is.list(theta)){
    theta<-theta
  }else if(theta=="start"){
    theta<-ergmm.fit$start
  }else if(theta=="mle"){
    theta<-summary(ergmm.fit,point.est=c("mle"),se=FALSE)$mle
  }else if(theta=="pmean"){
    theta<-summary(ergmm.fit,point.est=c("pmean"))$pmean
  }else if(theta=="mkl"){
    theta<-summary(ergmm.fit,point.est=c("mkl"))$mkl
  }else if(theta=="pmode"){
    theta<-summary(ergmm.fit,point.est=c("pmode"))$pmode
  }else if(is.numeric(theta) && round(theta)==theta){
    theta<-ergmm.fit$samples[[theta]]
  }else stop("Invalid parameter structure.")

  ergmm.EY.L(ergmm.fit$model,theta)
}
