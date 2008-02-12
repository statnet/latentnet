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



get.sample.deltas<-function(model,samples,control){
  ## Construct the "extended" beta: not just the coefficients, but also the scale and the average
  ## value of each random effect.
  beta.ext<-cbind(if(model$p) samples$beta, # covariate coefs
                  if(model$d) log(apply(sqrt(apply(samples$Z^2,1:2,sum)),1,mean)), # scale of Z
                  if(model$sender) apply(samples$sender,1,mean), # sender eff.
                  if(model$receiver) apply(samples$receiver,1,mean), # receiver eff.
                  if(model$sociality) apply(samples$sociality,1,mean)) # sociality eff.
  cov.beta.ext<-cov(beta.ext)
  
  
  ## Zero-out those rows in the variance-covariance matrix that are to be proposed independently.
  cov.pos<-0
  if(model$p){
    if("beta" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+(1:model$p))
    cov.pos<-cov.pos+model$p
  }
  
  if(model$d){
    if("Z.scl" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+1)
      cov.pos<-cov.pos+1
  }
  
  if(model$sender){
    if("sender" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+1)
    cov.pos<-cov.pos+1
  }
  
  if(model$receiver){
    if("receiver" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+1)
    cov.pos<-cov.pos+1
  }
  
  if(model$sociality){
    if("sociality" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+1)
    cov.pos<-cov.pos+1
  }
  
  ## Take the Choletsky decomposition of the empirical covariance matrix.
  control$group.deltas<-chol(cov.beta.ext)*control$pilot.factor
  control
}

## Zero the off-diagonal entries in a particular set of rows and columns.
partial.diag<-function(m,idx){
  m<-as.matrix(m)
  saved<-diag(m)
  m[idx,]<-0
  m[,idx]<-0
  diag(m)<-saved
  m
}
