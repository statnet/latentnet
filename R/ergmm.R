ergmm <- function(formula,response=NULL,family="Bernoulli.logit",fam.par=NULL,
                  control=ergmm.control(),
                  user.start=ergmm.par.blank(),
                  prior=as.ergmm.par.list(list(adjust.beta.var=TRUE)),
                  tofit=c("pmode","mcmc","mkl","mkl.mbc","mle","procrustes","klswitch"),
                  Z.ref=NULL,
                  Z.K.ref=NULL,
                  seed=NULL,
                  orthogonalize=FALSE,
                  verbose=FALSE){
  #    current.warn <- options()$warn
  #    options(warn=0)
  
  control$verbose<-verbose
  control$tofit<-ergmm.tofit.resolve(tofit)
  
  if(control$threads>1) require(snowFT)
  
  
  ## If the random seed has been specified, save the old seed, to
  ## pick up where it left off. If not, don't.
  if(!is.null(seed)){
    old.seed<-.Random.seed
    .Random.seed<-seed
  }else runif(1) # This is needed to initialize .Random.seed if it isn't already.
  start.seed<-.Random.seed
  
  if(class(formula)=="ergmm.model"){
    if(length(prior)){
      model<-formula
      prior<-prior
    }
    stop("If an ergmm.model is specified in place of a formula, prior must also be specified")
  }else{
    tmp <- ergmm.get.model(formula, response, family, fam.par, orthogonalize, prior)
    model<-tmp$model
    prior<-tmp$prior
  }
  burnin.start<-ergmm.initvals(model,user.start,prior,control)
  if(control$beta.delta.adjust) control<-adjust.beta.delta(model,control)
  
  if(control$tofit$mcmc){
    if(control$burnin>0){
      if(control$tune){
        if(control$verbose) cat("Tuning parameters for burnin...\n")
        tuning<-ergmm.tuner(model,burnin.start,prior,control)
        if(control$verbose) cat("Finished.\n")
        for(name in names(tuning)){
          control[[name]]<-tuning[[name]]
        }
      }

      if(control$verbose) cat("Burning in... ")
      
      if(control$threads<=1){
        # Burn in one thread.
        if(control$store.burnin){
          burnin.samples<-ergmm.MCMC.C(model,burnin.start,prior,control,
                                       samplesize=control$burnin/control$interval)$samples
          sampling.start<-burnin.samples[[length(burnin.samples$llk)]]
        }else sampling.start<-ergmm.MCMC.C(model,burnin.start,prior,control,
                                           samplesize=1,interval=control$burnin)$samples[[1]]
      }else{
        burnin.samples<-ergmm.MCMC.snowFT(control$threads,control$threads,
                                          model.l=list(model),
                                          start.l=list(burnin.start),
                                          prior.l=list(prior),
                                          control.l=list(control),
                                          samplesize.l=list(control$burnin/control$interval))$samples
        sampling.start<-sapply(1:control$threads,
                               function(thread) burnin.samples[[thread]][[length(burnin.samples[[thread]]$llk)]],
                               simplify=FALSE)
      }
      
      if(control$verbose) cat("Finished.\n")
    }else sampling.start<-burnin.start
    
    if(control$tune){
      if(control$verbose) cat("Tuning parameters for sampling run...\n ")
      tuning<-ergmm.tuner(model,sampling.start,prior,control,control$burnin>0 && control$threads>1)
      if(control$verbose) cat("Finished.\n")
      for(name in names(tuning)){
        control[[name]]<-tuning[[name]]
      }
    }
    
    if(control$verbose) cat("Starting sampling run... ")
    if(control$threads<=1)
      mcmc.out <-  ergmm.MCMC.C(model,sampling.start,prior,control)
    else{
      mcmc.out <- ergmm.MCMC.snowFT(control$threads,if(control$burnin) 1 else control$threads,
                                    model.l=list(model),
                                    start.l=if(control$burnin) sampling.start else list(sampling.start),
                                    prior.l=list(prior),
                                    control.l=list(control),
                                    samplesize.l=list(ceiling(control$samplesize/control$threads)))
      mcmc.out$samples <- stack.ergmm.par.list.list(mcmc.out$samples)
    }
    if(control$verbose) cat("Finished.\n")
  }
  else mcmc.out<-NULL
  
  v<-ergmm.statseval(mcmc.out, model, burnin.start,  prior, control,
                     Z.ref, Z.K.ref)
  
    if(control$tofit$mcmc){
      v$burnin.start<-burnin.start
      v$sampling.start<-sampling.start
      if(control$store.burnin)
        v$burnin.samples<-burnin.samples
    }
  
  v$starting.seed<-start.seed
  if(!is.null(seed)) .Random.seed<-old.seed
  
  v
}
