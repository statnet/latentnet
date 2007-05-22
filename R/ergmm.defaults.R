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
                        RE.delta=0.3,
                        RE.shift.delta=0.4,
                        beta.delta=0.4,
                        flyapart.penalty=0,
                        skip.MCMC=FALSE,
                        skip.Procrustes=FALSE,
                        skip.MBC=FALSE,
                        skip.KLswitch=FALSE,
                        store.burnin=FALSE,
                        verbose=FALSE){
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
       RE.delta=RE.delta,
       RE.shift.delta=RE.shift.delta,
       beta.delta=beta.delta,
       flyapart.penalty=flyapart.penalty,
       skip.MCMC=skip.MCMC,
       skip.Procrustes=skip.Procrustes,
       skip.MBC=skip.MBC,
       skip.KLswitch=skip.KLswitch,
       store.burnin=store.burnin,
       verbose=verbose)
}

