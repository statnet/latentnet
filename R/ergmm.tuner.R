get.init.deltas<-function(model, control){
  nterms<-model$p+(if(model$d) 1 else 0)+model$sender+model$receiver+model$sociality

  ## If proposal coefficient matrix is given, leave it alone.
  if(is.matrix(control$group.deltas)){
    if(any(dim(control$group.deltas)!=nterms))
      stop(paste("Incorrect proposal coefficient matrix size: (",
                 paste(dim(control$group.deltas),collapse=", "),"), ",
                 "while the model calls for ",nterms,".",sep=""))
    return(control)
  }
  
  ## If a vector of appropriate length is given, use a diagonal matrix
  if(length(control$group.deltas)==nterms){
    control$group.deltas<-diag(control$group.deltas)
    return(control)
  }
  
  ## If a scalar is given, construct a diagonal matrix that's in the ballpark.
  if(length(control$group.deltas)==1){
    control$group.deltas<-1/sapply(1:model$p,function(i) sqrt(mean((model$X[[i]][observed.dyads(model$Yg)])^2)))
    if(model$d) control$group.deltas<-c(control$group.deltas, 0.05)
    control$group.deltas<-c(control$group.deltas,rep(1,model$sender+model$receiver+model$sociality))
    control$group.deltas<-diag(control$group.deltas*control$group.deltas*2/(1+nterms),nrow=nterms)
  }

  control
}

get.sample.deltas<-function(model,sample,control){
  use.draws<-ceiling(length(sample)*control$pilot.discard.first):length(sample)
  if(model$d) control$Z.delta<-control$Z.delta*mean(sample$Z.rate[use.draws])/control$target.acc.rate
  if(model$sender || model$receiver || model$sociality) control$RE.delta<-control$RE.delta*mean(sample$Z.rate[use.draws])/control$target.acc.rate
  
  control$pilot.factor<-control$pilot.factor*mean(sample$beta.rate[use.draws])/control$target.acc.rate
  cov.beta.ext<-cov.beta.ext(model,sample[use.draws])
  
  ## Take the Choletsky decomposition of the empirical covariance matrix.
  control$group.deltas<-try(chol(cov.beta.ext)*control$pilot.factor)
  if(inherits(control$group.deltas,"try-error")) stop("Looks like a pilot run did not mix at all (practically no proposals were accepted). Try using a smaller starting proposal variance.")
  control
}

## Compute the empirical covariance of coefficients, latent scale, and random effect means.
cov.beta.ext<-function(model,sample){
  ## Construct the "extended" beta: not just the coefficients, but also the scale and the average
  ## value of each random effect.
  beta.ext<-cbind(if(model$p) sample$beta, # covariate coefs
                  if(model$d) log(apply(sqrt(apply(sample$Z^2,1:2,sum)),1,mean)), # scale of Z
                  if(model$sender) apply(sample$sender,1,mean), # sender eff.
                  if(model$receiver) apply(sample$receiver,1,mean), # receiver eff.
                  if(model$sociality) apply(sample$sociality,1,mean)) # sociality eff.
  cov(beta.ext)
}
