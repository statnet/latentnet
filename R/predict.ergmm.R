predict.ergmm<-function(object,...,type="post"){
  if(class(type)=="ergmm.par"){
    type<-type
  }else if(type=="start"){
    type<-object$start
  }else if(type=="mle"){
    type<-object$mle
  }else if(type=="pmean"){
    type<-summary(object,point.est=c("pmean"))$pmean
  }else if(type=="mkl"){
    type<-object$mkl
  }else if(type=="pmode"){
    type<-object$pmode
  }else if(type=="post"){
    return(with(object,post.predict.C(model,sample,control)))
  }else if(is.numeric(type) && round(type)==type){
    type<-object$sample[[type]]
  }else stop("Invalid parameter structure.")

  ergmm.EY(object$model,type,NA.unobserved=FALSE)
}
