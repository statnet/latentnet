ergmm <- function(formula,response=NULL,family="Bernoulli",fam.par=NULL,
                  control=ergmm.control(),
                  user.start=NULL,
                  prior=list(),
                  skipMCMC=FALSE,
                  skipProcrustes=FALSE,
                  skipMBC=FALSE,
                  Z.ref=NULL,
                  Z.K.ref=NULL,
                  store.burnin=FALSE,
                  randseed=NULL,
                  verbose=FALSE,
                  orthogonalize=FALSE)
{
#    current.warn <- options()$warn
#    options(warn=0)

    if(control$threads>1) require(snowFT)
    if(!is.null(randseed))
      set.seed(as.integer(randseed))

    tmp <- ergmm.build.model(formula, response, family, fam.par, orthogonalize, prior)
    model<-tmp$model
    prior<-tmp$prior

    mle1<-state<-ergmm.initvals(model,user.start,prior,control,verbose)

    if(!skipMCMC){
      if(control$burnin>0){
        if(control$tune){
          if(verbose) cat("Tuning parameters for burnin...\n")
          tuning<-ergmm.tuner(model,state,prior,control,verbose)
          if(verbose) cat("Finished.\n")
          for(name in names(tuning)){
            control[[name]]<-tuning[[name]]
          }
        }

        if(verbose) cat("Burning in... ")
        
        if(control$threads<=1){
          # Burn in one thread.
          if(store.burnin){
            burnin.samples<-ergmm.MCMC.C(model,state,prior,control,
                                         samplesize=control$burnin/control$interval)$samples
            state<-burnin.samples[[length(burnin.samples$llk)]]
          }
          else
            state<-ergmm.MCMC.C(model,state,prior,control,
                                samplesize=1,interval=control$burnin)$samples[[1]]
        }
        else{
          burnin.samples<-ergmm.MCMC.snowFT(control$threads,control$threads,
                                            model.l=list(model),
                                            start.l=list(state),
                                            prior.l=list(prior),
                                            control.l=list(control),
                                            samplesize.l=list(control$burnin/control$interval))$samples
          state<-sapply(1:control$threads,
                        function(thread) burnin.samples[[thread]][[length(burnin.samples[[thread]]$llk)]],
                        simplify=FALSE)
        }
        
        if(verbose) cat("Finished.\n")
      }
      else
        state<-mle1

      if(control$tune){
        if(verbose) cat("Tuning parameters for sampling run...\n ")
        tuning<-ergmm.tuner(model,state,prior,control,verbose,control$burnin>0 && control$threads>1)
        if(verbose) cat("Finished.\n")
        for(name in names(tuning)){
          control[[name]]<-tuning[[name]]
        }
      }
      
      if(verbose) cat("Starting sampling run... ")
      if(control$threads<=1)
        mcmc.out <-  ergmm.MCMC.C(model,state,prior,control)
      else{
        mcmc.out <- ergmm.MCMC.snowFT(control$threads,if(control$burnin) 1 else control$threads,
                                      model.l=list(model),
                                      start.l=if(control$burnin) state else list(state),
                                      prior.l=list(prior),
                                      control.l=list(control),
                                      samplesize.l=list(ceiling(control$samplesize/control$threads)))
        mcmc.out$samples <- stack.ergmm.par.list.list(mcmc.out$samples)
      }
      if(verbose) cat("Finished.\n")
    }
    else mcmc.out<-NULL
    
    v<-ergmm.statseval(mcmc.out, model, mle1,  prior, control,
                       skipProcrustes, Z.ref, Z.K.ref, verbose, skipMBC)

    if(!skipMCMC) {
      v$burnin.start<-mle1
      v$main.start<-state
      if(store.burnin)
        v$burnin.samples<-burnin.samples
    }
    
#    options(warn=current.warn)
    v
}
