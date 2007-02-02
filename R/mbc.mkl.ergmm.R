mbc.mkl.ergmm<-function(ergmm.fit,...){
  with(ergmm.fit,bayesmbc(model$G,mkl$Z,prior,...))$pmean
}

add.mkl.clusters.ergmm<-function(ergmm.fit,...){
  ergmm.fit$mkl<-c(ergmm.fit$mkl,mbc.mkl.ergmm(ergmm.fit,...))
  ergmm.fit
}
