predict.ergmm<-function(object,...,which.par="post"){
  if(class(which.par)=="ergmm.par"){
    which.par<-which.par
  }else if(which.par=="start"){
    which.par<-ergmm.fit$start
  }else if(which.par=="mle"){
    which.par<-ergmm.fit$mle
  }else if(which.par=="pmean"){
    which.par<-summary(ergmm.fit,point.est=c("pmean"))$pmean
  }else if(which.par=="mkl"){
    which.par<-ergmm.fit$mkl
  }else if(which.par=="pmode"){
    which.par<-ergmm.fit$pmode
  }else if(which.par=="post"){
    return(with(ergmm.fit,post.predict.C(model,samples,control)))
  }else if(is.numeric(which.par) && round(which.par)==which.par){
    which.par<-ergmm.fit$samples[[which.par]]
  }else stop("Invalid parameter structure.")

  ergmm.EY(ergmm.fit$model,which.par)
}
