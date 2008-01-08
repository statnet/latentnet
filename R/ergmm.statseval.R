ergmm.statseval <- function (mcmc.out, model, start, prior, control,Z.ref=NULL,Z.K.ref=NULL){
  if(control$verbose) cat("Post-processing the MCMC output:\n")
   
  Yg<- model$Yg
  G <- model$G
  d <- model$d
  n <- network.size(Yg)
  p <- model$p
  family <- model$family
  samplesize <- control$samplesize

  if(!is.null(Z.ref)){
    Z.ref<-scale(Z.ref,scale=FALSE)
    control$Z.ref<-Z.ref
  }
  if(!is.null(Z.K.ref)){
    control$Z.K.ref<-Z.K.ref
  }
  
  ## Start putting together the output structure.
  if(is.null(mcmc.out)) mcmc.out<-list()
  mcmc.out$model<-model
  mcmc.out$prior<-prior
  mcmc.out$control<-control
  mcmc.out$start<-start
  class(mcmc.out)<-"ergmm"


  if(control$tofit$klswitch) mcmc.out<-labelswitch.samples.ergmm(mcmc.out)
  if(control$tofit$pmode) mcmc.out<-add.mcmc.pmode.pmode.ergmm(mcmc.out)
  if(control$tofit$mle) mcmc.out<-add.mcmc.mle.mle.ergmm(mcmc.out)
  if(control$tofit$mkl) mcmc.out<-add.mkl.pos.ergmm(mcmc.out)
  if(control$tofit$mkl.mbc) mcmc.out<-add.mkl.mbc.ergmm(mcmc.out)
  if(control$tofit$procrustes) mcmc.out<-proc.samples.ergmm(mcmc.out)

  class(mcmc.out)<-"ergmm"
  return(mcmc.out)
}

statsreeval.ergmm<-function(x,Z.ref=NULL,Z.K.ref=NULL,rerun=FALSE){
  if(!is.null(Z.ref)){
    Z.ref<-scale(Z.ref,scale=FALSE)
    x$control$Z.ref<-Z.ref
  }
  if(!is.null(Z.K.ref)){
    x$control$Z.K.ref<-Z.K.ref
  }
  x<-labelswitch.samples.ergmm(x)
  if(rerun) x<-add.mcmc.pmode.pmode.ergmm(x)
  if(rerun) x<-add.mcmc.mle.mle.ergmm(x)
  if(rerun) x<-add.mkl.pos.ergmm(x)
  x<-add.mkl.mbc.ergmm(x)
  x<-proc.samples.ergmm(x)
  x
}

find.mkl<-function(model,samples,control){
  if(control$verbose>1) cat("Evaluating matrix of predicted dyad values and finding initial value... ")
  EY<-post.predict.C(model,samples,control,TRUE)
  EY[!observed.dyads(model$Yg)]<-NA
  if(control$verbose>1) cat("Finished.\nMaximizing...")

  mkl<-samples[[attr(EY,"s.MKL")]]
  model$Ym<-EY
  for(i in 1:control$mle.maxit){
    if(control$verbose>1) cat(i,"")
    mkl.old<-mkl
    mkl<-find.mle(model,start=mkl,control=control,mllk=FALSE)
    if(is.null(mkl)) stop("MKL failed!")
    if(isTRUE(all.equal(mkl.old,mkl))) break
  }
  mkl
}

nullapply<-function(X,margin,FUN,...){
  if(is.null(X)) return(X)
  else apply(X,margin,FUN,...)
}

add.mkl.mbc.ergmm<-function(x,Z.K.ref=best.avail.Z.K.ref.ergmm(x)){
  if(is.null(x$mkl) || x$model$d<=0 || x$model$G<=0){
    if(x$control$verbose) cat("MKL is not available or non-latent-cluster model. MKL MBC will not be fitted.\n")
    return(x)
  }else Z<-x$mkl$Z

  if(x$control$verbose) cat("Fitting MBC conditional on MKL locations... ")
  x$mkl$mbc<-bayesmbc(x$model$G,Z,x$prior,Z.K.ref,verbose=x$control$verbose)$pmean
  if(x$control$verbose) cat("Finished.\n")
  x
}

add.mcmc.mle.mle.ergmm<-function(x,Z.ref=best.avail.Z.ref.ergmm(x)){
  if(!is.null(x$start)){
    if(x$control$verbose) cat("Using the conditional posterior mode to seed an MLE fit... ")
    mle1<-find.mle.loop(x$model,x$start,control=x$control)
    if(x$control$verbose) cat(" Finished.\n")
  }else mle1<-list(llk=-Inf)
  
  if(!is.null(x$mcmc.mle)){
    ## Use the iteration with the highest probability to seed another
    ## shot at the MLE.
    
    if(x$control$verbose) cat("Using the highest-likelihood iteration to seed another MLE fit... ")
    mle2 <- find.mle.loop(x$model,x$mcmc.mle,control=x$control)
    if(x$model$d>0)
      mle2$Z<-scale(mle2$Z,scale=FALSE)
    if(x$control$verbose) cat("Finished.\n")
  }
  else mle2<-list(llk=-Inf)
  if(is.null(mle2)) mle2<-list(llk=-Inf)
  
  
  if(mle2$llk>mle1$llk) x$mle<-mle2 else x$mle<-mle1
  
  if(x$model$d>0){
    x$mle$Z<-scale(x$mle$Z,scale=FALSE)
    
    if(!require(shapes,quietly=TRUE)){
      stop("You need the 'shapes' package to summarize the fit of latent cluster models.")
    }
    x$mle$Z<-procOPA(Z.ref,x$mle$Z,scale=FALSE,reflect=TRUE)$Bhat
  }
  
  x
}

find.mle.loop<-function(model,start,control){
  mle<-start
  for(i in 1:control$mle.maxit){
    mle.old<-mle
    mle<-find.mle(model,mle,control=control)
    if(all.equal(mle.old,mle)[1]==TRUE) break
  }
  mle
}


add.mcmc.pmode.pmode.ergmm<-function(x,Z.ref=best.avail.Z.ref.ergmm(x)){
  if(!is.null(x$start)){
    if(x$control$verbose) cat("Double-checking conditional posterior mode estimate... ")
    pmode2<-find.pmode.loop(x$model,x$start,prior=x$prior,control=x$control)
    if(x$control$verbose) cat("Finished.\n")
  }else pmode2<-list(mlp=-Inf)
  
  if(!is.null(x$mcmc.pmode)){
    if(x$control$verbose) cat("Using MCMC posterior mode to seed another conditional posterior mode fit... ")
    x$pmode<-find.pmode.loop(x$model,x$mcmc.pmode,prior=x$prior,control=x$control)
    if(x$control$verbose) cat("Finished.\n")
    
  }
  if(pmode2$mlp>x$pmode$mlp) x$pmode<-pmode2
    
  if(x$model$d>0 && !is.null(x$pmode)){
    Z<-scale(x$pmode$Z,scale=FALSE)
    if(!require(shapes,quietly=TRUE)){
      stop("You need the 'shapes' package to summarize the fit of latent cluster models.")
    }
    P<-procOPA(Z.ref,Z,scale=FALSE,reflect=TRUE)$R
    x$pmode$Z<-x$pmode$Z%*%P
    if(!is.null(x$pmode$Z.mean))
      x$pmode$Z.mean<-x$pmode$Z.mean%*%P
  }
  x
}

find.pmode.loop<-function(model,start,prior,control){
  pmode<-start
  for(i in 1:control$mle.maxit){
    if(control$verbose>1) cat(i,"")
    pmode.old<-pmode
    pmode<-find.mpe(model,pmode,prior=prior,given=as.ergmm.par.list(list(Z.K=pmode$Z.K)),control=control)
    if(model$G>1) pmode$Z.K<-find.clusters(model$G,pmode$Z)$Z.K
    if(all.equal(pmode.old,pmode)[1]==TRUE) break
  }
  pmode
}

add.mkl.pos.ergmm<-function(x, Z.ref=best.avail.Z.ref.ergmm(x)){
  if(!is.null(x$samples) && x$model$d>0){
    if(x$control$verbose) cat("Fitting the MKL locations... ")
    x$mkl<-find.mkl(x$model,x$samples,x$control)
    if(!require(shapes,quietly=TRUE)){
      stop("You need the 'shapes' package to summarize the fit of latent cluster models.")
    }
  }
  if(!is.null(x$mkl$Z)) x$mkl$Z<-scale(x$mkl$Z,scale=FALSE)
  if(!is.null(x$mkl$Z)) x$mkl$Z<-procOPA(Z.ref,x$mkl$Z,scale=FALSE,reflect=TRUE)$Bhat
  if(x$control$verbose) cat("Finished.\n")
  x
}

proc.samples.ergmm<-function(x,Z.ref=best.avail.Z.ref.ergmm(x)){
  if(!is.null(x$samples) && x$model$d>0){
    if(x$control$verbose) cat("Performing Procrustes tranformation... ")
    x$samples<-proc.Z.mean.C(x$samples,Z.ref,verbose=x$control$verbose)
    if(x$control$verbose) cat("Finished.\n")
  }
  x
}

labelswitch.samples.ergmm<-function(x,Z.K.ref=best.avail.Z.K.ref.ergmm(x)){
  if(!is.null(x$samples) && x$model$G>1){
    if(x$control$verbose) cat("Performing label-switching... ")
    require(mclust,quiet=TRUE)
    Q.start<-switch.Q.K(Z.K.ref,x$model$G)
    x$samples<-klswitch.C(Q.start,x$samples,verbose=x$control$verbose)
    if(x$control$verbose) cat("Finished.\n")
  }
  x
}

best.avail.Z.ref.ergmm<-function(x){
  if(x$model$d==0) return(NULL)
  
  if(!is.null(x$control$Z.ref)) return(x$control$Z.ref)
  if(!is.null(x$mkl$Z)) return(x$mkl$Z)
  if(!is.null(x$pmode$Z)) return(x$pmode$Z)
  if(!is.null(x$mcmc.pmode$Z)) return(x$mcmc.pmode$Z)
  if(!is.null(x$mle$Z)) return(x$mle$Z)
  if(!is.null(x$mcmc.mle$Z)) return(x$mcmc.mle$Z)
  if(!is.null(x$start$Z)) return(x$start$Z)
}

best.avail.Z.K.ref.ergmm<-function(x){
  if(x$model$G==0) return(NULL)
  
  if(!is.null(x$control$Z.K.ref)) return(x$control$Z.K.ref)
  if(!is.null(attr(x$samples,"Q"))) return(apply(attr(x$samples,"Q"),1,which.max))
  if(!is.null(x$pmode$Z.K)) return(x$pmode$Z.K)
  if(!is.null(x$mkl$mbc$Z.K)) return(x$mkl$mbc$Z.K)
  if(!is.null(x$mcmc.pmode$Z.K)) return(x$mcmc.pmode$Z.K)
  if(!is.null(x$start$Z.K)) return(x$start$Z.K)
  
  return(find.clusters(x$model$G,best.avail.Z.ref.ergmm(x))$Z.K)
}
