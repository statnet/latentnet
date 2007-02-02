bayesmbc<-function(G,Z,prior,samplesize=5000,interval=10,burnin=1000){
  start<-find.clusters(G,Z)
  state<-start<-with(start,list(Z=Z,
                                Z.mean=Z.mean,
                                Z.var=Z.var,
                                Z.K=Z.K,
                                Z.pK=Z.pK)
                     )

  state<-bayesmbc.MCMC.C(G,start,prior,
                         samplesize=1,interval=burnin)$samples[[1]]
  state$Z<-Z
  mcmc.out<-bayesmbc.MCMC.C(G,state,prior,samplesize,interval)

  mcmc.out$samples<-do.label.switching(mcmc.out$samples,prior,start$Z.K,Z=Z)

  mcmc.out$pmean<-with(mcmc.out$samples,
                       list(Z.mean=apply(Z.mean,2:3,mean),
                            Z.var=apply(Z.var,2,mean),
                            Z.K=apply(Z.K,2,function(x)which.max(tabulate(x,G))),
                            Z.pZK=t(apply(Z.K,2,function(x)tabulate(x,G)/length(x)))
                            )
                       )
  mcmc.out
}
