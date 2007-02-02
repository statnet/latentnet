ergmm.statseval <- function (mcmc.out, model, start, prior, control,
                             skipProcrustes=FALSE, Z.ref=NULL, Z.K.ref=NULL, verbose=FALSE,
                             skipMBC=FALSE){
  Yg<- model$Yg
  G <- model$G
  d <- model$d
  n <- network.size(Yg)
  p <- model$p
  family <- model$family
  samplesize <- control$samplesize
  if(!is.null(mcmc.out)){
    ## Use the iteration with the highest probability to seed another
    ## shot at the MLE.
    
    if(verbose) cat("Using the highest-likelihood iteration to seed another MLE fit.\n")
    mle2 <- find.mle.L(model,mcmc.out$mcmc.mle,control=control,flyapart.penalty=control$flyapart.penalty)
    if(d>0)
      mle2$Z<-scale(mle2$Z,scale=FALSE)
  }
  else mle2<-list(llk=-Inf)
  if(is.null(mle2)) mle2<-list(llk=-Inf)
  
  start$llk<-ergmm.loglike.L(model,start)
  if(mle2$llk<start$llk){
    mle2$llk<-start$llk
    mle2$beta<-start$beta
    mle2$Z<-start$Z
    mle2$sender<-start$sender
    mle2$receiver<-start$receiver
    mle2$sociality<-start$sociality
  }
  if(!is.null(Z.ref)){
    Z.ref<-scale(Z.ref,scale=FALSE)
    if(!require(shapes,quietly=TRUE)){
      stop("You need the 'shapes' package to summarize the fit of latent cluster models.")
    }
    mle2$Z<-procOPA(Z.ref,mle2$Z,scale=FALSE,reflect=TRUE)$Bhat
  }

  if(!is.null(mcmc.out) && d>0){
    if(verbose) cat("Fitting the MKL locations... ")
    mkl<-find.mkl.L(model,mcmc.out,control)
    mkl$Z<-scale(mkl$Z,scale=FALSE)
    if(!is.null(Z.ref)){
      if(!require(shapes,quietly=TRUE)){
        stop("You need the 'shapes' package to summarize the fit of latent cluster models.")
      }
      mcmc.out$Z.ref<-Z.ref
      mkl$Z<-procOPA(Z.ref,mkl$Z,scale=FALSE,reflect=TRUE)$Bhat
    }
    if(verbose) cat("Finished.\n")
    if(G>0 && !skipMBC){
      if(verbose) cat("Fitting MKL clusters... ")
      mkl$mbc<-bayesmbc(G,mkl$Z,prior)$pmean
      if(!is.null(Z.K.ref)){
        perm<-nearest.perm(Z.K.ref,mkl$mbc$Z.K)
        mkl$mbc$Z.K<-ergmm.labelswitch(mkl$mbc$Z.K,perm)
        mkl$mbc$Z.mean<-mkl$mbc$Z.mu[perm,]
        mkl$mbc$Z.var<-mkl$mbc$Z.var[perm]
      }
      if(verbose) cat("Finished.\n")
    }
    mcmc.out$mkl<-mkl
  }
  
  if(!is.null(mcmc.out) && !skipProcrustes){
    ## Procrustify and label-switch.
    if(d>0){
      if(verbose) cat("Performing Procrustes tranformation... ")
      mcmc.out$samples<-proc.Z.mean(mcmc.out$samples,mkl$Z,center=FALSE)
      if(verbose) cat("Finished.\n")
    }
    if(G>1){
      clust2<-find.clusters(model$G,mkl$Z)
      if(!is.null(Z.K.ref)) {
        mcmc.out$Z.K.ref<-Z.K.ref
        clust2<-nearest.perm(Z.K.ref,clust2$Z.K)
      }
      else{
        clust2<-clust2$Z.K
      }
      if(verbose) cat("Performing label-switching... ")
      #mcmc.out$samples<-do.label.switching(mcmc.out$samples,mkl$Z,clust2)
      mcmc.out$samples<-do.label.switching(mcmc.out$samples,prior,clust2)
      if(verbose) cat("Finished.\n")
    }
  }

  ## Start putting together the output structure.
  mcmc.out$model<-model
  mcmc.out$prior<-prior
  mcmc.out$control<-control
  mcmc.out$mle<-mle2
  if(verbose) cat("Fitting another posterior mode estimate... ")
  mcmc.out$pmode<-find.mpe.L(model,mcmc.out$mcmc.pmode,prior=prior,control=control,flyapart.penalty=control$flyapart.penalty)
  cat("Finished.\b")
  
  class(mcmc.out)<-"ergmm"
  return(mcmc.out)
}

proc.Z.mean<-function(samples,Z.ref,center=FALSE){
  ## For procOPA
  if(!require(shapes,quietly=TRUE)){
    stop("You need the 'shapes' package to summarize the fit of latent cluster models.")
  }

  n<-dim(Z.ref)[1]
  G<-dim(samples$Z.mean)[2]
  k<-dim(Z.ref)[2]
  ## Center Z.ref.
  Z.ref<-Z.ref-cbind(rep(1,n))%*%nullapply(Z.ref,2,mean)

  for(i in 1:dim(samples$Z)[1]){

    if(center){
      Z.centroid<-nullapply(seldrop(samples$Z[i,,,drop=FALSE],1),2,mean)
      if(!is.null(samples$Z.mean))
        samples$Z.mean[i,,]<-seldrop(samples$Z.mean[i,,,drop=FALSE],1)-cbind(rep(1,G))%*%Z.centroid
      samples$Z[i,,]<-seldrop(samples$Z[i,,,drop=FALSE],1)-cbind(rep(1,n))%*%Z.centroid
    }
    
    Z<-seldrop(samples$Z[i,,,drop=FALSE],1)-cbind(rep(1,n))%*%nullapply(seldrop(samples$Z[i,,,drop=FALSE],1),2,mean)
    P<-procOPA(Z.ref,Z,FALSE,TRUE)$R
    samples$Z[i,,]<-seldrop(samples$Z[i,,,drop=FALSE],1)%*%P
    if(!is.null(samples$Z.mean))
      samples$Z.mean[i,,]<-seldrop(samples$Z.mean[i,,,drop=FALSE],1)%*%P
  }
  samples
}


do.label.switching<-function(samples,prior,Z.K.ref,Z=NULL){
  d<-dim(samples$Z.mean)[3]
  G<-dim(samples$Z.mean)[2]
  S<-dim(samples$Z.mean)[1]
  Z.n.ref<-tabulate(Z.K.ref)
  if(G>1){
    permute <- ergmm.permutation(G)
    for(s in 1:S){
      sample<-samples[[s]]
      if(is.null(Z)) Z<-samples$Z[s,,]
      
      Z.mean<-samples$Z.mean[s,,]
      Z.var<-samples$Z.var[s,]
      Z.K<-samples$Z.K[s,]
      Z.K.n<-c(tabulate(Z.K),numeric(G-max(Z.K)))
      
      Z.bar.ref<-nullapply(Z,2,function(x,Z.K.ref)tapply(x,Z.K.ref,mean),Z.K.ref = Z.K.ref)
      Z.bar<-Z.mean*(1+Z.var/prior$Z.mean.var*Z.K.n)
      
      best.perm<-permute[which.min(apply(permute,1,
                                         function(p){
                                           sum((Z.bar[p,]-Z.bar.ref)^2)
                                         })),]
      
      samples$Z.K[s,] <-ergmm.labelswitch(Z.K,best.perm)
      samples$Z.mean[s,,] <- Z.mean[best.perm,]
      samples$Z.var[s,] <- Z.var[best.perm]
    }
  }
  samples
}

do.label.switching.old<-function(samples,Z.ref,Z.K.ref){
  d<-dim(samples$Z.mean)[3]
  G<-dim(samples$Z.mean)[2]
  samplesize<-dim(samples$Z.mean)[1]
  if(G>1){
    mu.ref <- nullapply(Z.ref,2,function(x,ki)tapply(x,ki,mean),ki = Z.K.ref)
    
    permute <- ergmm.permutation(G)
    
    for(loop in 1:samplesize){
      sample<-samples[[loop]]
      mu.1 <- nullapply(if(!is.null(sample$Z))sample$Z else Z.ref,2,
                        function(x,ki)tapply(x,ki,mean),ki = sample$Z.K)
      mutab <- table(sample$Z.K)
      mu.names <- as.numeric(names(mutab))
      n1 <- length(mutab)
      d1 <- as.matrix(dist(rbind(mu.1,mu.ref)))[1:n1,(n1+1):(n1+G)]
      if(length(mutab)==1) {
        d1.min <- order(d1)[1]
        samples$Z.K[loop,] <- d1.min
        samples$Z.mean[loop,d1.min,] <- mu.1
      }
      else{
        d1.use <- rep(TRUE,G)
        if(length(mutab)<G){
          d1.use <- rep(FALSE,G)
          d1.new <- matrix(0,G,G)
          mu.new <- matrix(0,G,d) #should be #of dimensions
          j <- 1
          for(i in 1:G)
            if(any(mu.names == i)){
              d1.new[i,] <- d1[j,]
              mu.new[i,] <- mu.1[j,]
              j <- j + 1
              d1.use[i] <- TRUE
            }
          d1 <- d1.new
          mu.1 <- mu.new
        }
        
        d1.vec <- rep(0,nrow(permute))
        for(j in 1:nrow(permute))
          for(i in 1:G)
            d1.vec[j] <- d1.vec[j] + d1.use[i] * d1[i,permute[j,i]]
        
        d1.min <- order(d1.vec)[1]
        per.to <- order(permute[d1.min,])
        samples$Z.K[loop,] <-ergmm.labelswitch(samples$Z.K[loop,],per.to)
        samples$Z.mean[loop,,] <- mu.1[per.to,]
        samples$Z.var[loop,] <- samples$Z.var[loop,per.to]
      }
    }
  }
  return(samples)
}

find.mkl.L<-function(model,mcmc.out,control){
  EY.f<-mk.EY.f(model$family)
  EYm<-matrix(0,network.size(model$Yg),network.size(model$Yg))
  for(i in 1:control$samplesize){
    state<-mcmc.out$samples[[i]]
    eta<-ergmm.eta(model$Yg,
                   state$beta,
                   model$X,
                   state$Z,
                   state$sender,
                   state$receiver,
                   state$sociality)
    EYm<-EYm+EY.f(eta,model$fam.par)
  }
  EYm<-EYm/control$samplesize

  min.mkl<-control$samplesize
  min.SSY<-sum((EYm-EY.f(eta,model$fam.par))^2)
  for(i in (control$samplesize-1):1){
    state<-mcmc.out$samples[[i]]
    eta<-ergmm.eta(model$Yg,
                   state$beta,
                   model$X,
                   state$Z,
                   state$sender,
                   state$receiver,
                   state$sociality)
    SSY<-sum((EYm-EY.f(eta,model$fam.par))^2)
    if(SSY<min.SSY){
      min.SSY<-SSY
      min.mkl<-i
    }
  }

  find.mle.L(model,start=mcmc.out$samples[[min.mkl]],control=control,mllk=FALSE,Ym=EYm)
}

nullapply<-function(X,margin,FUN,...){
  if(is.null(X)) return(X)
  else apply(X,margin,FUN,...)
}

reeval.ergmm<-function(x){
  y<-ergmm.statseval(list(samples=x$samples,mcmc.mle=x$mcmc.mle),x$model,x$start,x$prior,x$control,Z.ref=NULL)
  if(!is.null(x$tuning)) y$tuning<-x$tuning
  if(!is.null(x$tuning.history)) y$tuning.history<-x$tuning.history
  if(!is.null(x$main.start)) y$main.start<-x$main.start
  y
}

relabel.ergmm<-function(x,per.to){
  for(s in 1:x$control$samplesize){
    if(!is.null(x$samples$Z.K))
      x$samples$Z.K[i,]<-ergmm.labelswitch(x$samples$Z.K[i,],per.to)
    if(!is.null(x$samples$Z.mean))
      x$samples$Z.mean[i,,]<-x$samples$Z.mean[i,per.to,]
    if(!is.null(x$samples$Z.var))
      x$samples$Z.var[i,]<-x$samples$Z.var[i,per.to]
  }
  invisible(x)
}

reProcrustify.ergmm<-function(x,Z.ref,center=FALSE){
  x$mcmc.mle$Z<-procOPA(Z.ref,x$mcmc.mle$Z,scale=FALSE,reflect=TRUE)$Bhat
  x$mkl$Z<-procOPA(Z.ref,x$mkl$Z,scale=FALSE,reflect=TRUE)$Bhat
  x$samples<-proc.Z.mean(x$samples,x$mkl$Z,center=center)
  invisible(x)
}
