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

mvimode <- function(x, fix=NULL, maxit=10, theta0=NULL, simple=TRUE)
{
 if(!is.matrix(x)){
  stop("mvimode() requires a sample matrix")
 }
 theta <- apply(x,2,median)
 if(nrow(x) < 100 | simple){
  return(list(theta=theta,loglikelihoodratio=0))
 }
 if(is.null(fix)){
   fix <- rep(FALSE,ncol(x))
 }
 if(!require(locfit,quietly=TRUE, warn.conflicts = FALSE)){
   stop("You need the 'locfit' package to fit latent models.")
 }
#tempout <- capture.output(
# require(locfit, quietly = TRUE, warn.conflicts = FALSE)
#                         )
 if(!require(mvtnorm,quietly=TRUE, warn.conflicts = FALSE)){
   stop("You need the 'mvtnorm' package to fit latent cluster models.")
 }
 ncolx <- ncol(x) - sum(fix)
 x <- x[,!fix]
 covar0 <- var(x)
 stheta <- sqrt(diag(covar0))
 ntheta <- theta[!fix]
 if(ncolx==1){
#  1 variable
   x.den <- locfit.raw(x=x, alpha=0.7, maxk=1000, scale=T)
   ntheta <- x[which.max(predict(x.den, newdata=x))]
   loglik <- log(predict(x.den, newdata=ntheta, what="coef"))
   if(!is.null(theta0)){
    loglik <- loglik - log(predict(x.den, newdata=theta0, what="coef"))
   }
  }
 if(ncolx == 2)
 {
   for(j in 1:maxit){
    for(i in 1:ncol(x)){
     distval <- (x[,-i] - ntheta[-i])/sqrt(covar0[-i,-i])
     logdet <- log(covar0[-i,-i])
     logretval <- -(log(2 * pi) + logdet + distval*distval)/2
     ccc <- exp(logretval)
     ddd <- rank(-ccc)[1:min(nrow(x),200)]
     x.den <- locfit.raw(x=x[ddd,i], alpha=0.9, maxk=1000, weights=ccc[ddd])
     ccc <- try(predict(x.den, newdata=x[ddd,i])) 
     if(!inherits(ccc,"try-error")){
      ntheta[i] <- x[ddd,i][which.max(ccc)] 
     }
    }
   }
   loglik <- log(predict(x.den, newdata=ntheta[i], what="coef"))
   if(!is.null(theta0)){
    loglik <- loglik - log(predict(x.den, newdata=theta0[i], what="coef"))
   }
 }
 if(ncolx > 2)
 {
   for(j in 1:maxit){
    for(i in 1:ncol(x)){
     distval <- try(mahalanobis(x[,-i], center = ntheta[-i],
                    cov = covar0[-i,-i], tol.inv=1e-12))
     if(!inherits(distval,"try-error")){
      logdet <- sum(log(eigen(as.matrix(covar0[-i,-i]), symmetric = TRUE, 
                       only.values = TRUE)$values))
      logretval <- -(ncol(x[,-i]) * log(2 * pi) + logdet + distval)/2
      ccc <- exp(logretval)
      ddd <- rank(-ccc)[1:min(nrow(x),200)]
      x.den <- locfit.raw(x=x[ddd,i], alpha=0.9, maxk=1000, weights=ccc[ddd])
      ccc <- try(predict(x.den, newdata=x[ddd,i])) 
      if(!inherits(ccc,"try-error")){
       ntheta[i] <- x[ddd,i][which.max(ccc)] 
      }
     }
    }
   }
   loglik <- log(predict(x.den, newdata=ntheta[i], what="coef"))
   if(!is.null(theta0)){
    loglik <- loglik - log(predict(x.den, newdata=theta0[i], what="coef"))
   }
 }
 theta[!fix] <- ntheta
 list(theta=theta,loglikelihoodratio=loglik)
}

# Retrieve the number of free dyads (i.e., number of non-missing) of network x.
#
network.dyadcount<-function(x){
  if(!is.network(x))
    stop("network.dyadcount requires an argument of class network.")

  nodes <- network.size(x)
  if(is.directed(x)){
     dyads <- nodes * (nodes-1)
  }else{
   if(is.bipartite(x)){
    nevent <- get.network.attribute(x,"bipartite")
    nactor <- nodes - nevent
    dyads <- nactor * nevent
   }else{
    dyads <- nodes * (nodes-1)/2
   }
  }
#
# Adjust for missing
#
  design <- get.network.attribute(x,"design")
  if(!is.null(design)){
   dyads <- dyads - network.edgecount(design)
  }
  dyads
}
