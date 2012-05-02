proc.Z.mean.C<-function(sample,Z.ref,center=FALSE,verbose=0){
  n<-dim(Z.ref)[1]
  G<-dim(sample[["Z.mean"]])[2]
  if(is.null(G)) G<-0
  d<-dim(Z.ref)[2]
  S<-dim(sample[["Z"]])[1]
  ## Center Z.ref.
  Z.ref<-scale(Z.ref,scale=FALSE)

  Cret<-.C("procr_transform_wrapper",
           S=as.integer(S),
           n=as.integer(n),
           d=as.integer(d),
           G=as.integer(G),
           Z.ref=as.double(Z.ref),
           Z=as.double(sample[["Z"]]),
           Z.mean=as.double(sample[["Z.mean"]]),
           verbose=as.integer(verbose),
           
           PACKAGE="latentnet")
  sample[["Z"]]<-if(d>0)array(Cret[["Z"]],dim=c(S,n,d))
  sample[["Z.mean"]]<-if(G>0)array(Cret[["Z.mean"]],dim=c(S,G,d))
  
  sample
}
