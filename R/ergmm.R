#  File R/ergmm.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

## If the random seed has been specified, save the old seed, to
## pick up where it left off. If not, don't.

.save_set_seed <- function(seed){
  if(!is.null(seed)){
    old.seed<-try(.Random.seed)
    if(inherits(old.seed,"try-error")) old.seed<-NULL
    if(length(seed)==1) set.seed(seed)
    else .Random.seed<<-as.integer(seed)
  }else{
    runif(1) # This is needed to initialize .Random.seed if it isn't already.
    old.seed <- NULL
  }

  old.seed
}

.restore_set_seed <- function(old.seed){
  if(!is.null(old.seed)) .Random.seed<<-old.seed
}

#' Fit a Latent Space Random Graph Model
#' 
#' \code{\link{ergmm}} is used to fit latent space and latent space cluster
#' random network models, as described by Hoff, Raftery and Handcock (2002),
#' Handcock, Raftery and Tantrum (2005), and Krivitsky, Handcock, Raftery, and
#' Hoff (2009).  \code{\link{ergmm}} can return either a Bayesian model fit or
#' the two-stage MLE.
#' 
#' 
#' @aliases ergmm latent latentcluster
#' @param formula An formula object, of the form \code{g ~ <term 1> + <term 2>
#' ...}, where \code{g} is a network object or a matrix that can be coerced to
#' a network object, and \code{<term 1>}, \code{<term 2>}, etc., are each terms
#' for the model. See \code{\link{terms.ergmm}} for the terms that can be
#' fitted.  To create a network object in , use the \code{network} function,
#' then add nodal attributes to it using \code{set.vertex.attribute} if
#' necessary.
#' 
#' Note that, as in \code{\link{lm}}, the model will include an
#' \code{\link{intercept}} term. This behavior can be overridden by including a
#' \code{-1} or \code{+0} term in the formula, and a
#' \code{\link[=terms.ergmm]{1(mean=...,var=...)}} term can be used to set a
#' prior different from default.
#' @param response An optional edge attribute that serves as the response
#' variable. By default, presence (1) or absence (0) of an edge in \code{g} is
#' used.
#' @param family A character vector specifying the conditional distribution of
#' each edge value. See \link{families.ergmm} for the currently implemented
#' families.
#' @param fam.par For those families that require additional parameters, a
#' list.
#' @param control The MCMC parameters that do not affect the posterior
#' distribution such as the sample size, the proposal variances, and tuning
#' parameters, in the form of a named list. See \code{\link{control.ergmm}} for
#' more information and defaults.
#' @param user.start An optional initial configuration parameters for MCMC in
#' the form of a list. By default, posterior mode conditioned on cluster
#' assignments is used. It is permitted to only supply some of the parameters
#' of a configuration. If this is done, the remaining paramters are fitted
#' conditional on those supplied.
#' @param prior The prior parameters for the model being fitted in the form of
#' a named list. See \link{terms.ergmm} for the names to use.  If given, will
#' override those given in the formula terms, making it useful as a convenient
#' way to store and reproduce a prior distribution. The list or prior
#' parameters can also be extracted from an \link[=ergmm.object]{ERGMM fit
#' object}. See \code{\link{ergmm.prior}} for more information.
#' @param tofit A character vector listing some subset of "pmode", "mcmc",
#' "mkl", "mkl.mbc", "mle","procrustes", and "klswitch", defaulting to all of
#' the above, instructing \code{\link{ergmm}} what should be returned as a part
#' of the \link[=ergmm.object]{ERGMM fit object}. Omiting can be used to skip
#' particular steps in the fitting process. If the requested procedure or
#' output depends on some other procedure or output not explictly listed, the
#' dependency will be resolved automatically.
#' @param Z.ref If given, used as a reference for Procrustes analysis.
#' @param Z.K.ref If given, used as a reference for label-switching.
#' @param seed If supplied, random number seed.
#' @param verbose If this is \code{TRUE} (or \code{1}), causes information to
#' be printed out about the progress of the fitting, particularly initial value
#' generation. Higher values lead to greater verbosity.
#' @return \code{\link{ergmm}} returns an object of class
#' \code{\link[=ergmm.object]{ergmm}} containing the information about the
#' posterior.
#' @seealso network, set.vertex.attributes, set.network.attributes,
#' summary.ergmm, terms.ergmm, families.ergmm
#' @references Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum (2002).
#' \emph{Model-Based Clustering for Social Networks.} Journal of the Royal
#' Statistical Society: Series A, 170(2), 301-354.
#' 
#' Peter D. Hoff, Adrian E. Raftery and Mark S. Handcock (2002).  \emph{Latent
#' space approaches to social network analysis.} Journal of the American
#' Statistical Association, 97(460), 1090-1098.
#' 
#' Pavel N. Krivitsky, Mark S. Handcock, Adrian E. Raftery, and Peter D. Hoff
#' (2009). \emph{Representing degree distributions, clustering, and homophily
#' in social networks with latent cluster random effects models}.  Social
#' Networks, 31(3), 204-213.
#' 
#' Pavel N. Krivitsky and Mark S. Handcock (2008).  \emph{Fitting Position
#' Latent Cluster Models for Social Networks with \code{latentnet}}. Journal of
#' Statistical Software, 24(5).
#' @keywords graphs
#' @examples
#' 
#' \donttest{
#' #
#' # Use 'data(package = "latentnet")' to list the data sets in a
#' #
#' data(package="latentnet")
#' #
#' # Using Sampson's Monk data, lets fit a 
#' # simple latent position model
#' #
#' data(sampson)
#' samp.fit <- ergmm(samplike ~ euclidean(d=2))
#' #
#' # See if we have convergence in the MCMC
#' mcmc.diagnostics(samp.fit)
#' #
#' # Plot the fit
#' #
#' plot(samp.fit)
#' #
#' # Using Sampson's Monk data, lets fit a latent clustering random effects model
#' #
#' samp.fit2 <- ergmm(samplike ~ euclidean(d=2, G=3)+rreceiver)
#' #
#' # See if we have convergence in the MCMC
#' mcmc.diagnostics(samp.fit2)
#' #
#' # Plot the fit.
#' #
#' plot(samp.fit2, pie=TRUE)
#' }
#' @import network ergm
#' @export ergmm
ergmm <- function(formula,response=NULL,family="Bernoulli",fam.par=NULL,
                  control=control.ergmm(),
                  user.start=list(),
                  prior=ergmm.prior(),
                  tofit=c("mcmc","mkl","mkl.mbc","procrustes","klswitch"),
                  Z.ref=NULL,
                  Z.K.ref=NULL,
                  seed=NULL,
                  verbose=FALSE){
  
  control[["verbose"]]<-verbose
  control[["tofit"]]<-ergmm.tofit.resolve(tofit)
  
  if(control[["threads"]]>1||control[["kl.threads"]]>1){
    with(control,
         if(sample.size%%threads || (burnin/interval)%%threads)
         stop("Please make the MCMC sample size and the ratio burnin/interval a multiple of the number of threads."))
    if(!requireNamespace("snowFT", quietly=TRUE))
      stop("Package 'snowFT' is required for multithreaded MCMC.")
  }

  old.seed <- .save_set_seed(seed)
  
  start.seed<-.Random.seed
  
  if(class(formula)=="ergmm.model"){
    if(length(prior)){
      model<-formula
      prior<-prior
    }
    stop("If an ergmm.model is specified in place of a formula, prior must also be specified")
  }else{
    tmp <- ergmm.get.model(formula, response, family, fam.par, prior)
    model<-tmp[["model"]]
    prior<-tmp[["prior"]]
  }

  burnin.start<-burnin.state<-ergmm.initvals(model,user.start,prior,control)
  if(control[["tofit"]][["mcmc"]]){
    burnin.control<-get.init.deltas(model, control)
    burnin.controls<-list()
    burnin.samples<-list()

    if(control[["burnin"]]>0){
      burnin.runs<-max(control[["pilot.runs"]],1)
      burnin.size<-burnin.control[["burnin"]]/burnin.runs/burnin.control[["interval"]]

      if(control[["threads"]]>1) burnin.state<-list(burnin.state)
      if(burnin.control[["verbose"]]) cat("Burning in... ")
      for(pilot.run in 1:burnin.runs){
        if(burnin.control[["verbose"]]>1) cat(pilot.run,"")
        ## Set up a loop such that if a pilot run is catastrophically
        ## bad (only accepts or rejects a very tiny fraction of
        ## proposals), the proposals are "backed off" and the pilot
        ## run is redone.
        backoff<-TRUE
        while(backoff){
          burnin.controls[[length(burnin.controls)+1]]<-burnin.control
        
          if(burnin.control[["threads"]]<=1){
            ## Burn in one thread.
            burnin.sample<-ergmm.MCMC.C(model,burnin.state,prior,burnin.control,
                                        sample.size=burnin.size)[["sample"]]
            burnin.state<-burnin.sample[[burnin.size]]
          }else{
            ## Burn in multiple threads.
            burnin.sample<-ergmm.MCMC.snowFT(burnin.control[["threads"]],burnin.control[["threads"]],
                                             model.l=list(model),
                                             start.l=burnin.state,
                                             prior.l=list(prior),
                                             control.l=list(burnin.control),
                                             sample.size.l=list(burnin.size))[["sample"]]
            burnin.state<-sapply(1:burnin.control[["threads"]],
                                 function(thread) burnin.sample[[thread]][[burnin.size]],
                                 simplify=FALSE)
            burnin.sample<-.stack.ergmm.par.list.list(burnin.sample)
          }
          if(control[["store.burnin"]]) burnin.samples[[length(burnin.samples)+1]]<-burnin.sample
          if(control[["pilot.runs"]]){
            burnin.control<-backoff.check(model,burnin.sample,burnin.control)
            backoff<-burnin.control[["backedoff"]]
          }else backoff<-FALSE
        }
        
        if(control[["pilot.runs"]]) burnin.control<-get.sample.deltas(model, burnin.sample, burnin.control)
      }
      if(burnin.control[["verbose"]]) cat("Finished.\n")
    }
    
    sampling.start<-burnin.state
    
    control<-burnin.control
    
    if(control[["verbose"]]) cat("Starting sampling run... ")
    if(control[["threads"]]<=1)
      mcmc.out <-  ergmm.MCMC.C(model,sampling.start,prior,control)
    else{
      mcmc.out <- ergmm.MCMC.snowFT(control[["threads"]],if(control[["burnin"]]) 1 else control[["threads"]],
                                    model.l=list(model),
                                    start.l=if(control[["burnin"]]) sampling.start else list(sampling.start),
                                    prior.l=list(prior),
                                    control.l=list(control),
                                    sample.size.l=list(ceiling(control[["sample.size"]]/control[["threads"]])))
      mcmc.out[["sample"]] <- .stack.ergmm.par.list.list(mcmc.out[["sample"]])
    }
    if(control[["verbose"]]) cat("Finished.\n")
  }
  else mcmc.out<-NULL
  
  v<-ergmm.statseval(mcmc.out, model, burnin.start,  prior, control,
                     Z.ref, Z.K.ref)
  
    if(control[["tofit"]][["mcmc"]]){
      if(control[["burnin"]]){
        v[["burnin.start"]]<-burnin.start
        v[["burnin.controls"]]<-burnin.controls
        if(control[["store.burnin"]]) v[["burnin.samples"]]<-burnin.samples
      }
      v[["sampling.start"]]<-sampling.start
    }
  
  v[["starting.seed"]]<-start.seed

  .restore_set_seed(old.seed)
  
  v
}
