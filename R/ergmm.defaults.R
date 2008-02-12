## These functions define defaults passed to ergmm(...).
## They are meant to be used in a manner similar to glm.control(...)
## for glm.

## Variables affecting the sampling process but not the posterior distribution.
ergmm.control<-function(samplesize=2000,
                        burnin=1000,
                        interval=10,
                        threads=1,
                        mle.maxit=400,
                        Z.delta=0.4,
                        RE.delta=0.3,
                        propose.ind=c(),
                        pilot.runs=2,
                        pilot.factor=0.3,
                        group.deltas=0.3,
                        accept.all=FALSE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}

ergmm.fit.deps<-list(pmode=character(0),
                     mcmc=character(0),
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
