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
