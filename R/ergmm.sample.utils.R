proc.Z.mean.C<-function(samples,Z.ref,center=FALSE,verbose=0){
  n<-dim(Z.ref)[1]
  G<-dim(samples$Z.mean)[2]
  d<-dim(Z.ref)[2]
  S<-dim(samples$Z)[1]
  ## Center Z.ref.
  Z.ref<-scale(Z.ref,scale=FALSE)

  Cret<-.C("procr_transform_wrapper",
           S=as.integer(S),
           n=as.integer(n),
           d=as.integer(d),
           G=as.integer(G),
           Z.ref=as.double(Z.ref),
           Z=as.double(samples$Z),
           Z.mean=as.double(samples$Z.mean),
           verbose=as.integer(verbose),
           PACKAGE="latentnet")
  samples$Z<-if(d>0)array(Cret$Z,dim=c(S,n,d))
  samples$Z.mean<-if(G>0)array(Cret$Z.mean,dim=c(S,G,d))
  
  samples
}

klswitch.C <- function(Q.start,samples,Z=NULL,maxit=100,verbose=0)
{
  
  Z.ref<-!is.null(Z)

  n <- if(Z.ref) dim(Z)[1] else dim(samples$Z)[2]
  
  d <- dim(samples$Z.mean)[3]
  S <- dim(samples$Z.mean)[1]
  G <- dim(samples$Z.mean)[2]

  if(!all(dim(Q.start)==c(n,G))) stop("Incorrect dimensions for initial Q matrix.")
  
  Cret <- .C("klswitch_wrapper",
             maxit = as.integer(maxit),
             S = as.integer(S),
             n = as.integer(n),
             d = as.integer(d),
             G = as.integer(G),
             Z = if(Z.ref) as.double(Z) else as.double(samples$Z),
             Z.ref=as.integer(Z.ref),
             Z.mean = as.double(samples$Z.mean),
             Z.var = as.double(samples$Z.var),
             Z.K = as.integer(samples$Z.K),
             Z.pK = as.double(samples$Z.pK),
             Q = as.double(Q.start),
             verbose=as.integer(verbose),
             PACKAGE="latentnet")
  
  samples$Z.mean<-array(Cret$Z.mean,dim=c(S,G,d))
  samples$Z.var<-matrix(Cret$Z.var,S,G)
  samples$Z.K<-matrix(Cret$Z.K,S,n)
  samples$Z.pK<-matrix(Cret$Z.K,S,G)
  attr(samples,"Q")<-array(Cret$Q,dim=c(1,n,G))[1,,]

  samples
}

thin.ergmm<-function(x,by){
  if(x$control$threads>1) warning("Multithreaded run output. Stuff might be broken.")
  S<-x$control$samplesize
  s.kept<-seq(from=1,to=S,by=by)
  x$samples<-x$samples[s.kept]
  x$control$interval<-x$control$interval*by
  x$control$samplesize<-length(s.kept)
  x
}

switch.Q.K<-function(K,G,smooth=1/G){
  n<-length(K)
  Q<-matrix(smooth,n,G)
  for(i in 1:n) Q[i,K[i]]<-1+smooth
  t(apply(Q,1,function(x) x/sum(x)))
}

post.predict.C<-function(model,samples,control){
  n<-network.size(model$Yg)
  d<-model$d
  p<-model$p

  ## Figure out the design matrix.
  observed<-observed.dyads(model$Yg)

  if((observed==(diag(n)==0) && is.directed(model$Yg)) ||
     (observed==lower.tri(diag(n)) && !is.directed(model$Yg)))
    observed<-NULL

  familyID<-switch(model$family,
                   Bernoulli=0,
                   binomial=1,
                   Poisson=2)
  
  array(.C("post_pred_wrapper",
           S = as.integer(control$samplesize),
           
           n = as.integer(n),
           p = as.integer(p),
           d = as.integer(d),
           
           dir=is.directed(model$Yg),
           family=as.integer(familyID),
           iconsts=as.integer(model$iconsts),
           dconsts=as.double(model$dconsts),
           
           X=as.double(unlist(model$X)),
           
           Z = as.double(samples$Z),
           beta = as.double(samples$beta), # coef
           sender = as.double(samples$sender),
           receiver = as.integer(samples$receiver),
           sociality = as.integer(model$sociality),
           observed=as.integer(observed),
           
           EY=double(n*n),
           
           verbose=as.integer(control$verbose),
           PACKAGE="latentnet")$EY,dim=c(1,n,n))[1,,] 
}
