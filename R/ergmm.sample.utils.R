proc.Z.mean.C<-function(sample,Z.ref,center=FALSE,verbose=0){
  n<-dim(Z.ref)[1]
  G<-dim(sample$Z.mean)[2]
  if(is.null(G)) G<-0
  d<-dim(Z.ref)[2]
  S<-dim(sample$Z)[1]
  ## Center Z.ref.
  Z.ref<-scale(Z.ref,scale=FALSE)

  Cret<-.C("procr_transform_wrapper",
           S=as.integer(S),
           n=as.integer(n),
           d=as.integer(d),
           G=as.integer(G),
           Z.ref=as.double(Z.ref),
           Z=as.double(sample$Z),
           Z.mean=as.double(sample$Z.mean),
           verbose=as.integer(verbose),
           PACKAGE="latentnet")
  sample$Z<-if(d>0)array(Cret$Z,dim=c(S,n,d))
  sample$Z.mean<-if(G>0)array(Cret$Z.mean,dim=c(S,G,d))
  
  sample
}

klswitch.C <- function(Q.start,sample,Z=NULL,maxit=100,verbose=0)
{
  
  Z.ref<-!is.null(Z)

  n <- if(Z.ref) dim(Z)[1] else dim(sample$Z)[2]
  
  d <- dim(sample$Z.mean)[3]
  S <- dim(sample$Z.mean)[1]
  G <- dim(sample$Z.mean)[2]

  if(!all(dim(Q.start)==c(n,G))) stop("Incorrect dimensions for initial Q matrix.")
  
  Cret <- .C("klswitch_wrapper",
             maxit = as.integer(maxit),
             S = as.integer(S),
             n = as.integer(n),
             d = as.integer(d),
             G = as.integer(G),
             Z = if(Z.ref) as.double(Z) else as.double(sample$Z),
             Z.ref=as.integer(Z.ref),
             Z.mean = as.double(sample$Z.mean),
             Z.var = as.double(sample$Z.var),
             Z.K = as.integer(sample$Z.K),
             Z.pK = as.double(sample$Z.pK),
             Q = as.double(Q.start),
             verbose=as.integer(verbose),
             PACKAGE="latentnet")
  
  sample$Z.mean<-array(Cret$Z.mean,dim=c(S,G,d))
  sample$Z.var<-matrix(Cret$Z.var,S,G)
  sample$Z.K<-matrix(Cret$Z.K,S,n)
  sample$Z.pK<-matrix(Cret$Z.pK,S,G)
  attr(sample,"Q")<-array(Cret$Q,dim=c(1,n,G))[1,,]

  sample
}

labelswitch.K <- function(Z.K,perm)
{
  Z.K.new <- Z.K
  for(i in 1:length(perm))
    Z.K.new[Z.K == perm[i]] <- i
  return(Z.K.new)
}

klswitch.snowFT<-function(threads,Q.start,sample,Z=NULL,maxit=100,verbose=0){
  Z.ref<-!is.null(Z)
  n <- if(Z.ref) dim(Z)[1] else dim(sample$Z)[2]
  d <- dim(sample$Z.mean)[3]
  S <- dim(sample$Z.mean)[1]
  G <- dim(sample$Z.mean)[2]

  if(!all(dim(Q.start)==c(n,G))) stop("Incorrect dimensions for initial Q matrix.")
  

  if(!require(snowFT)) stop("Package snowFT required for multithreaded KL Switching!")
  
  Cret <- .C("klswitch_pK_wrapper",
             S = as.integer(S),
             n = as.integer(n),
             d = as.integer(d),
             G = as.integer(G),
             Z = if(Z.ref) as.double(Z) else as.double(sample$Z),
             Z.ref = as.integer(Z.ref),
             Z.mean = as.double(sample$Z.mean),
             Z.var = as.double(sample$Z.var),
             Z.K = as.integer(sample$Z.K),
             Z.pK = as.double(sample$Z.pK),
             verbose = as.integer(verbose),
             pK = double(S*n*G),
             PACKAGE="latentnet")

  pK<-array(Cret$pK,dim=c(S,n,G))

  pK.l<-list()
  for(thread in 1:threads){
    pK.l[[thread]]<-pK[(1+(thread-1)*S/threads):(thread*S/threads),,]
  }

  Q<-Q.start
  
  if(verbose>0) cat("KLswitch: Iterating between label-switching to Q and recalculating Q.\n")
  for(it in 1:maxit){
    changed<-FALSE

    best.perms.l<-{
      if(threads==1)
        list(klswitch.step2.snowFT.slave(1,lib=path.to.me,Q=Q,pK.l=pK.l))
      else performParallel(threads,1:threads,
                           klswitch.step2.snowFT.slave,
                           lib=path.to.me,
                           Q=Q,
                           pK.l=pK.l)
    }
    
    Z.K<-sample$Z.K
    Z.pK<-sample$Z.pK
    Z.mean<-sample$Z.mean
    Z.var<-sample$Z.var
    
    for(thread in 1:threads){
      from.s<-1+(thread-1)*S/threads
      to.s<-thread*S/threads
      pK<-pK.l[[thread]]
      for(s in from.s:to.s){
        best.perm<-best.perms.l[[thread]][s-from.s+1,]
        if(all(best.perm==0)) next

        changed<-TRUE

        pK[s-from.s+1,,]<-pK[s-from.s+1,,best.perm]
        Z.pK[s,]<-Z.pK[s,best.perm]
        Z.mean[s,,]<-Z.mean[s,best.perm,]
        Z.var[s,]<-Z.var[s,best.perm]
        Z.K[s,]<-labelswitch.K(Z.K[s,],best.perm)
      }
      pK.l[[thread]]<-pK
    }
    
    sample$Z.K<-Z.K
    sample$Z.pK<-Z.pK
    sample$Z.mean<-Z.mean
    sample$Z.var<-Z.var

    Q<-t(apply(Z.K,2,tabulate,nbins=G))/S
    if(verbose>2)
      cat("KLswitch: Iterating: Completed ",it,"/",maxit,".\n",sep="")
    if(!changed){
      if(verbose>1) cat("KLswitch: Converged.\n")
      break
    }
  }

  attr(sample,"Q")<-Q
  
  sample
}

klswitch.step2.snowFT.slave<-function(i,lib,Q,pK.l){
  library(latentnet,lib=lib)

  n<-dim(Q)[1]
  G<-dim(Q)[2]
  S<-dim(pK.l[[i]])[1]
  
  Cret <- .C("klswitch_step2_wrapper",
             S = as.integer(S),
             n = as.integer(n),
             G = as.integer(G),
             Q = as.double(Q),
             pK = as.double(pK.l[[i]]),
             best.perms = integer(S*G),
             PACKAGE="latentnet")

  matrix(Cret$best.perms,nrow=S,ncol=G)
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

switch.Q.K<-function(K,G,smooth=1/G){
  n<-length(K)
  Q<-matrix(smooth,n,G)
  for(i in 1:n) Q[i,K[i]]<-1+smooth
  t(apply(Q,1,function(x) x/sum(x)))
}

post.predict.C<-function(model,sample,control,MKL=FALSE){
  n<-network.size(model$Yg)
  d<-model$d
  p<-model$p

  ## Figure out the design matrix.
  observed<-observed.dyads(model$Yg)

  if((observed==(diag(n)==0) && is.directed(model$Yg)) ||
     (observed==lower.tri(diag(n)) && !is.directed(model$Yg)))
    observed<-NULL

  ret<-.C("post_pred_wrapper",
          S = as.integer(control$sample.size),
          
          n = as.integer(n),
          p = as.integer(p),
          d = as.integer(d),
          
          dir=is.directed(model$Yg),
          family=as.integer(model$familyID),
          iconsts=as.integer(model$iconsts),
          dconsts=as.double(model$dconsts),
          
          X=as.double(unlist(model$X)),
          
          Z = as.double(sample$Z),
          beta = as.double(sample$beta), # coef
          sender = as.double(sample$sender),
          receiver = as.double(sample$receiver),
          sociality = as.double(model$sociality),
          observed=as.integer(observed),
          
          EY=double(n*n),
          s.MKL=if(MKL) integer(1) else integer(0),
          verbose=as.integer(control$verbose),
          PACKAGE="latentnet")
  EY<-array(ret$EY,dim=c(1,n,n))[1,,] 
  if(MKL) attr(EY,"s.MKL")<-ret$s.MKL+1 # C counts from 0; R counts from 1
  EY
}

post.predict.R<-function(model,sample,control,MKL=FALSE){
  EY.f<-EY.fs[[model$familyID]]
  EY<-matrix(0,network.size(model$Yg),network.size(model$Yg))
  for(i in 1:control$sample.size){
    state<-sample[[i]]
    eta<-ergmm.eta(model,state)
    EY<-EY+EY.f(eta,model$fam.par)
  }
  EY<-EY/control$sample.size

  if(MKL){
    min.MKL<-NA
    min.dev<-Inf
    model$Ym<-EY
    model$Ym[!observed.dyads(model$Yg)]<-NA
    for(i in 1:control$sample.size){
      state<-sample[[i]]
      dev<--ergmm.lpY(model,state,up.to.const=TRUE)
      if(dev<min.dev){
        min.dev<-dev
        min.MKL<-i
      }
    }
    attr(EY,"s.MKL")<-min.MKL
  }
  
  EY
}
