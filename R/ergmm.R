ergmm <- function(formula,response=NULL,family="Bernoulli.logit",fam.par=NULL,
                  control=ergmm.control(),
                  user.start=ergmm.par.blank(),
                  prior=ergmm.prior(),
                  tofit=c("mcmc","mkl","mkl.mbc","procrustes","klswitch"),
                  Z.ref=NULL,
                  Z.K.ref=NULL,
                  seed=NULL,
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
    tmp <- ergmm.get.model(formula, response, family, fam.par, prior)
    model<-tmp$model
    prior<-tmp$prior
  }
  
  if(control$tofit$mcmc){
    burnin.start<-burnin.state<-ergmm.initvals(model,user.start,prior,control)
    burnin.control<-get.init.deltas(model, control)
    burnin.controls<-list()

    if(control$burnin>0){
      burnin.runs<-max(control$pilot.runs,1)
      burnin.size<-burnin.control$burnin/burnin.runs/burnin.control$interval

      if(control$threads>1) burnin.state<-list(burnin.state)
      if(burnin.control$verbose) cat("Burning in... ")
      for(pilot.run in 1:control$pilot.runs){
        if(burnin.control$verbose>1) cat(pilot.run,"")
        burnin.controls[[length(burnin.controls)+1]]<-burnin.control
        
        if(burnin.control$threads<=1){
          ## Burn in one thread.
          burnin.samples<-ergmm.MCMC.C(model,burnin.state,prior,burnin.control,
                                       samplesize=burnin.size)$samples
          burnin.state<-burnin.samples[[burnin.size]]
        }else{
          ## Burn in multiple threads.
          burnin.samples<-ergmm.MCMC.snowFT(burnin.control$threads,burnin.control$threads,
                                            model.l=list(model),
                                            start.l=burnin.state,
                                            prior.l=list(prior),
                                            control.l=list(burnin.control),
                                            samplesize.l=list(burnin.size))$samples
          burnin.state<-sapply(1:burnin.control$threads,
                               function(thread) burnin.samples[[thread]][[burnin.size]],
                               simplify=FALSE)
          burnin.samples<-stack.ergmm.par.list.list(burnin.samples)

        }
        if(control$pilot.runs) burnin.control<-get.sample.deltas(model, burnin.samples, burnin.control)
      }
      if(burnin.control$verbose) cat("Finished.\n")
    }
    
    sampling.start<-burnin.state
    
    control<-burnin.control
    
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
      if(control$burnin){
        v$burnin.start<-burnin.start
        v$burnin.controls<-burnin.controls
        v$burnin.samples<-burnin.samples
      }
      v$sampling.start<-sampling.start
    }
  
  v$starting.seed<-start.seed
  if(!is.null(seed)) .Random.seed<-old.seed
  
  v
}
