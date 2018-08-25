#  File R/ergmm.defaults.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
## These functions define defaults passed to ergmm(...).
## They are meant to be used in a manner similar to glm.control(...)
## for glm.

## Variables affecting the sampling process but not the posterior distribution.


#' Auxiliary for Controlling ERGMM Fitting
#' 
#' Auxiliary function as user interface for \code{ergmm} fitting. Typically
#' only used when calling \code{ergmm}. It is used to set parameters that
#' affect the sampling but do not affect the posterior distribution.
#' 
#' 
#' @param sample.size The number of draws to be taken from the posterior
#' distribution.
#' @param burnin The number of initial MCMC iterations to be discarded.
#' @param interval The number of iterations between consecutive draws.
#' @param threads The number of chains to run. If greater than 1, package
#' \code{\link[snowFT:snowFT-package]{snowFT}} is used to take advantage of any
#' multiprocessing or distributed computing capabilities that may be available.
#' Currently, only PVM (via \code{rpvm}) has been tested. Note, also, that PVM
#' daemon needs to be started before the package is loaded.
#' @param kl.threads If greather than 1, uses an experimental parallelized
#' label-switching algorithm. This is not guaranteed to work.
#' @param mle.maxit Maximum number of iterations for computing the starting
#' values, posterior modes, MLEs, MKL estimates, etc..
#' @param Z.delta Standard deviation of the proposal for the jump in the
#' individual latent space position, or its starting value for the tuner.
#' @param RE.delta Standard deviation of the proposal for the jump in the
#' individual random effects values, or its starting value for the tuner.
#' @param group.deltas A scalar, a vector, or a matrix of an appropriate size,
#' giving the initial proposal structure for the ``group proposal'' of a jump
#' in covariate coefficients, scaling of latent space positions, and a shift in
#' random ffects. If a matrix of an appropriate size is given, it is used as a
#' matrix of coefficients for a correlated proposal. If a vector is given, an
#' independent proposal is used with the corresponding elements being proposal
#' standard deviations. If a scalar is given, it is used as a multiplier for an
#' initial heuristic for the proposal structure. It is usually best to leave
#' this argument alone and let the adaptive sampling be used.
#' @param pilot.runs Number of pilot runs into which to split the burn-in
#' period. After each pilot run, the proposal standard deviations and
#' coefficients \code{Z.delta}, \code{RE.delta}, and \code{group.deltas} are
#' reevaluated. If set to \code{0}, disables adaptive sampling, and only makes
#' a single burn-in run.
#' @param pilot.factor Initial value for the factor by which the coefficients
#' gotten by a Choletsky decomposition of the pilot sample covariance matrix
#' are multiplied.
#' @param pilot.discard.first Proportion of draws from the pilot run to discard
#' for estimating acceptance rate and group proposal covariance.
#' @param target.acc.rate Taget acceptance rate for the proposals used. After a
#' pilot run, the proposal variances are adjusted upward if the acceptance rate
#' is above this, and downward if below.
#' @param backoff.threshold If a pilot run's acceptance rate is below this,
#' redo it with drastically reduced proposal standard deviation.  Set to
#' \code{0} to disable this behavior.
#' @param backoff.factor Factor by which to multiply the relevant proposal
#' standard deviation if its acceptance rate fell below the backoff threshold.
#' @param accept.all Forces all proposals to be accepted unconditionally. Use
#' only in debugging proposal distributions!
#' @param store.burnin If \code{TRUE}, the samples from the burnin are also
#' stored and returned, to be used in MCMC diagnostics.
#' @param refine.user.start If \code{TRUE}, the values passed to
#' \code{\link{ergmm}} in the \code{user.start} argument can be updated by the
#' mode-finding algorithm.
#' @return A list with the arguments as components.
#' @seealso \code{\link{ergmm}}
#' @keywords graphs misc
#' @examples
#' 
#' \donttest{
#' data(sampson)
#' ## Shorter run than default.
#' samp.fit<-ergmm(samplike~euclidean(d=2,G=3)+rreceiver,
#' control=ergmm.control(burnin=1000,sample.size= 2000,interval=5))
#' }
#' 
#' @export
control.ergmm<-function(sample.size=4000,
                        burnin=10000,
                        interval=10,
                        threads=1,
                        kl.threads=1,
                        mle.maxit=100,
                        Z.delta=0.6,
                        RE.delta=0.6,
                        group.deltas=0.4,
                        pilot.runs=4,
                        pilot.factor=0.8,
                        pilot.discard.first=0.5,
                        target.acc.rate=0.234,
                        backoff.threshold=0.05,
                        backoff.factor=0.2,
                        accept.all=FALSE,
                        store.burnin=FALSE,
                        refine.user.start=TRUE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}

#' @rdname control.ergmm
#' @export
ergmm.control <- control.ergmm

#' Auxiliary for Setting the ERGMM Prior
#' 
#' Auxiliary function as user interface for \code{\link{ergmm}} prior
#' specification. Typically only used when calling \code{\link{ergmm}}. It is
#' used to supply the parameters of the prior distribution of the model, to
#' overwrite those specified in the model formula, and to supply miscellaneous
#' prior parameters.
#' 
#' 
#' @param \dots Prior distribution parameters. See \link{terms.ergmm} for more
#' information.
#' @param adjust.beta.var A shortcut: whether the prior variance for each
#' covariate coefficient should be divided by the mean square of that
#' covariate. This adjustment affects those variances specified in the formula
#' or by default, but not those specified through the \code{prior=} argument.
#' @return A list with the arguments as elements.
#' @seealso \code{\link{ergmm}}, \code{\link{terms.ergmm}}
#' @keywords graphs models
#' @export ergmm.prior
ergmm.prior<-function(...,adjust.beta.var=TRUE){
  prior<-list(...)
  prior[["adjust.beta.var"]]<-adjust.beta.var
  as.ergmm.par.list(prior)
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
