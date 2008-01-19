## These functions define defaults passed to ergmm(...).
## They are meant to be used in a manner similar to glm.control(...)
## for glm.

## Variables affecting the sampling process but not the posterior distribution.
ergmm.control<-function(samplesize=2000,
                        burnin=1000,
                        interval=10,
                        threads=1,
                        mle.maxit=400,
                        tune=FALSE,
                        tuning.runs=100,
                        tuning.runsize=8,
                        Z.delta=0.4,
                        Z.tr.delta=0.4,
                        Z.scl.delta=0.02,
                        beta.delta=0.4,
                        store.burnin=FALSE){
  list(samplesize=samplesize,
       burnin=burnin,interval=interval,
       threads=threads,
       mle.maxit=mle.maxit,
       tune=tune,
       tuning.runs=tuning.runs,
       tuning.runsize=tuning.runsize,
       Z.delta=Z.delta,
       Z.tr.delta=Z.tr.delta,
       Z.scl.delta=Z.scl.delta,
       beta.delta=beta.delta,
       store.burnin=store.burnin)
}

ergmm.fit.deps<-list(pmode=character(0),
                     mcmc=c("pmode"),
                     mkl=c("mcmc"),
                     mkl.mbc=c("mkl"),
                     mle=c("pmode"),
                     klswitch=c("mcmc"),
                     procrustes=c("mcmc"))

ergmm.tofit.resolve<-function(tofit){
  if(class(tofit)=="list"){
    tofit.c<-c()
    for(fit in names(tofit))
      if(tofit[[fit]])
        tofit.c<-c(tofit.c,fit)
    tofit<-tofit.c
  }
  
  oldlen<-length(tofit)-1
  while(length(tofit)!=oldlen){
    oldlen<-length(tofit)
    for(fit in tofit)
      tofit<-c(tofit,ergmm.fit.deps[[fit]])
    tofit<-unique(tolower(tofit))
  }
  
  tofit.l<-list()
  for(fit in names(ergmm.fit.deps))
    tofit.l[[fit]] <- fit %in% tofit
  
  tofit.l
}
