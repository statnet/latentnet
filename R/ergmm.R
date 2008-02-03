ergmm <- function(formula,response=NULL,family="Bernoulli.logit",fam.par=NULL,
                  control=ergmm.control(),
                  user.start=ergmm.par.blank(),
                  prior=as.ergmm.par.list(list(adjust.beta.var=TRUE)),
                  tofit=c("pmode","mcmc","mkl","mkl.mbc","mle","procrustes","klswitch"),
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
  burnin.start<-ergmm.initvals(model,user.start,prior,control)
   
  if(control$tofit$mcmc){
    
    if(control$burnin>0){
      burnin.control<-get.group.deltas(control$group.deltas, model, NULL, control)

      if(burnin.control$verbose) cat("Burning in... ")
      if(burnin.control$threads<=1){
        # Burn in one thread.
        burnin.samples<-ergmm.MCMC.C(model,burnin.start,prior,burnin.control,
                                     samplesize=burnin.control$burnin/burnin.control$interval)$samples
        sampling.start<-burnin.samples[[length(burnin.samples$llk)]]
      }else{
        burnin.samples<-ergmm.MCMC.snowFT(burnin.control$threads,burnin.control$threads,
                                          model.l=list(model),
                                          start.l=list(burnin.start),
                                          prior.l=list(prior),
                                          control.l=list(burnin.control),
                                          samplesize.l=list(burnin.control$burnin/burnin.control$interval))$samples
        sampling.start<-sapply(1:burnin.control$threads,
                               function(thread) burnin.samples[[thread]][[length(burnin.samples[[thread]]$llk)]],
                               simplify=FALSE)
        burnin.samples<-stack.ergmm.par.list.list(burnin.samples)
      }
      if(burnin.control$verbose) cat("Finished.\n")
    }else sampling.start<-burnin.start

    control<-get.group.deltas(control$group.deltas, model, if(control$pilot.runs && control$burnin) burnin.samples, control)
    
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
        v$burnin.control<-burnin.control
        v$burnin.samples<-burnin.samples
      }
      v$sampling.start<-sampling.start
    }
  
  v$starting.seed<-start.seed
  if(!is.null(seed)) .Random.seed<-old.seed
  
  v
}
