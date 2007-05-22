bayesmbc<-function(G,Z,prior,Z.K.ref=NULL,samplesize=2000,interval=10,burnin=500,verbose=FALSE){
  start<-find.clusters(G,Z)
  state<-start<-with(start,list(Z=Z,
                                Z.mean=Z.mean,
                                Z.var=Z.var,
                                Z.K=Z.K,
                                Z.pK=Z.pK)
                     )
  if(verbose>1) cat("Running MBC MCMC... ")
  state<-bayesmbc.MCMC.C(G,start,prior,
                         samplesize=1,interval=burnin)$samples[[1]]
  state$Z<-Z
  mcmc.out<-bayesmbc.MCMC.C(G,state,prior,samplesize,interval)
  if(verbose>1) cat("Finished.\n")

  if(verbose>1) cat("Running label switching... ")
  Q.start <- t(apply(unmap(if(is.null(Z.K.ref)) start$Z.K else Z.K.ref)+1/G,1,function(x) x/sum(x)))
  mcmc.out$samples<-klswitch.C(Q.start,mcmc.out$samples,Z)
  if(verbose>1) cat("Finished.\n")
  
  mcmc.out$pmean<-with(mcmc.out$samples,
                       list(Z.mean=apply(Z.mean,2:3,mean),
                            Z.var=apply(Z.var,2,mean),
                            Z.K=apply(Z.K,2,function(x)which.max(tabulate(x,G))),
                            Z.pZK=t(apply(Z.K,2,function(x)tabulate(x,G)/length(x)))
                            )
                       )
  mcmc.out
}
