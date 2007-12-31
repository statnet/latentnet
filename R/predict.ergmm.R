predict.ergmm<-function(object,...,which.par="post"){
  if(class(which.par)=="ergmm.par"){
    which.par<-which.par
  }else if(which.par=="start"){
    which.par<-object$start
  }else if(which.par=="mle"){
    which.par<-object$mle
  }else if(which.par=="pmean"){
    which.par<-summary(object,point.est=c("pmean"))$pmean
  }else if(which.par=="mkl"){
    which.par<-object$mkl
  }else if(which.par=="pmode"){
    which.par<-object$pmode
  }else if(which.par=="post"){
    return(with(object,post.predict.C(model,samples,control)))
  }else if(is.numeric(which.par) && round(which.par)==which.par){
    which.par<-object$samples[[which.par]]
  }else stop("Invalid parameter structure.")

  ergmm.EY(object$model,which.par)
}
