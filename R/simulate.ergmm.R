#  File R/simulate.ergmm.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
#' Draw from the distribution of an Exponential Random Graph Mixed Model
#' 
#' If passed a \code{\link[=ergmm.object]{ergmm}} fit object, \code{simulate}
#' is used to simulate networks from the posterior of an exponetial random
#' graph mixed model fit. Alternatively, a
#' \code{ergmm.model} can be passed to
#' simulate based on a particular parametr configuration.
#' 
#' A sample of networks is randomly drawn from the specified model. If a needed
#' value of \code{par} is missing, it is generated from its prior distribution.
#' 
#' @aliases simulate simulate.ergmm.model simulate.ergmm
#' @param object either an object of class \code{\link[=ergmm.object]{ergmm}}
#' for posterior simulation, or an object of class
#' \code{ergmm.model} for a specific model.
#' @param nsim number of networks to draw (independently)
#' @param seed random seed to use; defaults to using the current state of the
#' random number generator
#' @param par a list with the parameter configuration based on which to
#' simulate
#' @param prior a list with the prior distribution parameters that deviate from
#' their defaults
#' @param \dots Additional arguments. Currently unused.
#' @return If \code{nsim = 1}, \code{simulate} returns an object of class
#' \code{\link[network]{network}}. Otherwise, an object of class
#' \code{network.series} that is a list consisting of the following elements:
#' \item{\$formula}{The formula used to generate the sample.}
#' \item{\$networks}{A list of the generated networks.}
#' @seealso \code{\link{ergmm}}, \code{ network},
#' \code{\link[network]{print.network}}
#' @keywords graphs models nonlinear nonparametric cluster datagen
#' @examples
#' 
#' #
#' # Fit a short MCMC run: just the MCMC.
#' #
#' data(sampson)
#' gest <- ergmm(samplike ~ euclidean(d=2,G=3),
#'               control=ergmm.control(burnin=100,interval=5,sample.size=100),tofit="mcmc")
#' #
#' # Draw from the posterior
#' #
#' g.sim <- simulate(gest)
#' plot(g.sim)
#' #
#' # Draw from the first draw from the posterior
#' #
#' g.sim <- with(gest,simulate(model,par=sample[[1]],prior=prior))
#' plot(g.sim)
#' @importFrom stats simulate
#' @export
simulate.ergmm<-function(object, nsim=1, seed=NULL,...){
  extraneous.argcheck(...)

  old.seed <- .save_set_seed(seed)
  
  l<-list()
  for(i in seq_len(nsim)){
    iter<-floor(runif(1,1,object[["control"]][["sample.size"]]+1))
    l[[i]]<-sim.1.ergmm(object[["model"]],object[["sample"]][[iter]],object[["prior"]])
  }
  
  .restore_set_seed(old.seed)

  if(nsim > 1){
    l <- list(formula = object[["model"]][["formula"]], networks = l,
                     stats = NULL, coef=NULL)
    attr(l,"class")<-"network.list"
  }else{
    l <- l[[1]]
  }
  return(l)
}

#' @rdname simulate.ergmm
#' @export
simulate.ergmm.model<-function(object,nsim=1,seed=NULL,par,prior=list(),...){
  extraneous.argcheck(...)

  old.seed <- .save_set_seed(seed)

  l<-list()
  for(i in seq_len(nsim)){
    l[[i]]<-sim.1.ergmm(object,par,prior)
  }
  
  .restore_set_seed(old.seed)

  if(nsim==1) return(l[[1]])
  else{
    attr(l,"class")<-"network.list"
    return(l)
  }
}

#' @importFrom stats rchisq
sim.1.ergmm<-function(model,par,prior=list()){
  nv<-network.size(model[["Yg"]])
  mypar<-par
  
  if(length(model[["X"]])>0 && is.null(mypar[["beta"]]))
    mypar[["beta"]]<-rnorm(length(model[["X"]]),prior[["beta.mean"]],sqrt(prior[["beta.var"]]))

  if(model[["d"]]>0 && is.null(mypar[["Z"]])){
    if(model[["G"]]>0){
      if(is.null(mypar[["Z.mean"]]))
        mypar[["Z.mean"]]<-matrix(rnorm(model[["G"]]*model[["d"]],0,sqrt(prior[["Z.mean.var"]])),nrow=model[["G"]])
      if(is.null(mypar[["Z.K"]]))
        mypar[["Z.K"]]<-sample(seq_len(model[["G"]]),nv,replace=TRUE)
    }
    
    if(is.null(mypar[["Z.var"]]))
      mypar[["Z.var"]]<-prior[["Z.var"]]*prior[["Z.var.df"]]/rchisq(max(model[["G"]],1),prior[["Z.var.df"]])

    mypar[["Z"]]<-matrix(rnorm(nv*model[["d"]],
                            if(model[["G"]]>0) mypar[["Z.mean"]][mypar[["Z.K"]],] else 0,
                            if(model[["G"]]>0) mypar[["Z.var"]][mypar[["Z.K"]]] else mypar[["Z.var"]]
                            ),nrow=nv)
  }

  if(model[["sociality"]] && is.null(mypar[["sociality"]])){
    if(is.null(mypar[["sociality.var"]]))
      mypar[["sociality.var"]]<-with(prior,sociality.var*sociality.var.df/rchisq(1,sociality.var.df))
    model[["sociality"]]<-rnorm(nv,0,sqrt(mypar[["sociality.var"]]))
  }

  if(model[["sender"]] && is.null(mypar[["sender"]])){
    if(is.null(mypar[["sender.var"]]))
      mypar[["sender.var"]]<-with(prior,sender.var*sender.var.df/rchisq(1,sender.var.df))
    mypar[["sender"]]<-rnorm(nv,0,sqrt(mypar[["sender.var"]]))
  }
  
  if(model[["receiver"]] && is.null(mypar[["receiver"]])){
    if(is.null(mypar[["receiver.var"]]))
      mypar[["receiver.var"]]<-with(prior,receiver.var*receiver.var.df/rchisq(1,receiver.var.df))
    mypar[["receiver"]]<-rnorm(nv,0,sqrt(mypar[["receiver.var"]]))
  }

  if(model[["dispersion"]] && is.null(mypar[["dispersion"]])){
    mypar[["dispersion"]]<-with(prior,dispersion*dispersion.df/rchisq(1,dispersion.df))
  }

  eta<-ergmm.eta(model,mypar)

  sm<-rsm.fs[[model[["familyID"]]]](eta,dispersion=mypar[["dispersion"]],fam.par=model[["fam.par"]])
  if(!has.loops(model[["Yg"]]))
    diag(sm)<-0
  net<-with(model,as.network.matrix(sm,matrix.type="adjacency",
                                    directed=is.directed(Yg),
                                    loops=has.loops(Yg)))
  if(!is.null(model[["response"]]))
    net<-set.edge.value(net,model[["response"]],sm)

  attr(net,"ergmm.par")<-mypar

  net
}
