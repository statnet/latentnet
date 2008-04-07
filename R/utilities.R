#
# Return TRUE iff object x is a ergmm fit object
# or a latent model
#
is.latent<-function(x) inherits(x,"ergmm") && x$model$d>0

#
# Return TRUE iff object x is a ergmm fit object
# or a latent model and a "latentcluster" model or fit
#
is.latent.cluster<-function(x) inherits(x,"ergmm") && x$model$d>0 && x$model$G>0

#
# Return TRUE iff object x is a ergmm fit object
# or a latent model and a sender random effect
#
is.sender<-function(x) inherits(x,"ergmm") && x$model$sender
is.receiver<-function(x) inherits(x,"ergmm") && x$model$receiver

extraneous.argcheck<-function(...){
  ## Because #$^$% R wants various functions implementing generics to
  ## functions to take ..., which is wonderful for missing spelling
  ## errors.
  
  if(length(list(...)))stop("Extraneous arguments passed: ",
                            paste(list(...)))
}


thin.ergmm<-function(x,by){
  if(x$control$threads>1) warning("Multithreaded run output. Stuff might be broken.")
  S<-x$control$sample.size
  s.kept<-seq(from=1,to=S,by=by)
  x$sample<-x$sample[s.kept]
  x$control$interval<-x$control$interval*by
  x$control$sample.size<-length(s.kept)
  x
}

xtabs.ergmm<-function(x,ref){
  ref->Reference
  apply(attr(x$sample,"Q"),1,which.max)->Fitted
  xtabs(~Reference+Fitted)
}
