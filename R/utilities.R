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

sociomatrix<-function(object, matrix.type="adjacency", attrname=NULL)
{
  if(class(object)=="network.series") {
    object <- object$networks[[1]]
  }
  if(is.network(object)){
    x <- as.matrix.network(object,matrix.type=matrix.type, attrname=attrname)
    xnames <- network.vertex.names(object)
    dimnames(x) <- list(xnames,xnames)
    if(is.bipartite(object)){
     nevents <- is.bipartite(object)
     nactors <- network.size(object) - nevents
     x <- x[1:nactors, (1:nevents)+nactors]
    }
    x
  }
}

extraneous.argcheck<-function(...){
  ## Because #$^$% R wants various functions implementing generics to
  ## functions to take ..., which is wonderful for missing spelling
  ## errors.
  
  if(length(list(...)))stop("Extraneous arguments passed: ",
                            paste(list(...)))
}
