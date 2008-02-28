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

extraneous.argcheck<-function(...){
  ## Because #$^$% R wants various functions implementing generics to
  ## functions to take ..., which is wonderful for missing spelling
  ## errors.
  
  if(length(list(...)))stop("Extraneous arguments passed: ",
                            paste(list(...)))
}
