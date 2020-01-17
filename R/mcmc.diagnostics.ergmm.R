#  File R/mcmc.diagnostics.ergmm.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
#' Conduct MCMC diagnostics on an ERGMM fit
#' 
#' This function creates simple diagnostic plots for the MCMC sampled
#' statistics produced from a fit. It also prints the Raftery-Lewis
#' diagnostics, indicates if they are sufficient, and suggests the run length
#' required.
#' 
#' Produces the plots per \code{which.diags}.  Autocorrelation function that is
#' printed if "acf" is requested is for lags \code{0} and \code{interval}.
#' 
#' @aliases mcmc.diagnostics.ergmm mcmc.diagnostics
#' @param object An object of class \code{\link[=ergmm.object]{ergmm}}.
#' @param which.diags A list of diagnostics to produce. "cor" is the
#' correlation matrix of the statistics, "acf" plots the autocorrelation
#' functions, "trace" produces trace plots and density estimates, and "raftery"
#' produces the Raftery-Lewis statistics.
#' @param burnin If not \code{FALSE}, rather than perform diagnostics on the
#' sampling run, performs them on the pilot run whose index is given.
#' @param which.vars A named list mapping variable names to the indices to
#' include. If given, overrides the defaults and all arguments that follow.
#' @param vertex.i A numeric vector of vertices whose latent space coordinates
#' and random effects to include.
#' @param \dots Additional arguments. None are supported at the moment.
#' @return \code{mcmc.diagnostics.ergmm} returns a table of Raftery-Lewis
#' diagnostics.
#' @seealso \code{\link{ergmm}}, \code{\link{ergmm.object}},
#' \code{\link[coda]{raftery.diag}}, \code{\link[coda]{autocorr}},
#' \code{\link[coda]{plot.mcmc.list}}
#' @keywords graphs hplot debugging
#' @examples
#' 
#' \donttest{
#' #
#' data(sampson)
#' #
#' # test the mcmc.diagnostics function
#' #
#' gest <- ergmm(samplike ~ euclidean(d=2),
#'               control=ergmm.control(burnin=1000,interval=5))
#' summary(gest)
#' #
#' # Plot the traces and densities
#' #
#' mcmc.diagnostics(gest)
#' }
#' 
#' @export
mcmc.diagnostics.ergmm <- function(object,which.diags=c("cor","acf","trace","raftery"),
                                   burnin=FALSE,
                                   which.vars=NULL,
                                   vertex.i=c(1),...){
  extraneous.argcheck(...)
  
  if(is.null(object[["sample"]])) stop("MCMC was not run for this ERGMM fit.")

  x <- as.mcmc.list.ergmm(object,burnin,which.vars,vertex.i)
  oldask=par("ask")
  on.exit(par(ask=oldask))
  #' @importFrom grDevices dev.interactive
  par(ask=dev.interactive())

  if("cor" %in% which.diags){
    #' @importFrom coda autocorr autocorr.plot
    x.ac<-autocorr(x,lags=0:1)
    for(chain in seq(along=x.ac)){
      cat(paste("Chain",chain,"\n"))
      didnt.mix<-colnames(x.ac[[chain]][2,,])[which(is.nan(diag(x.ac[[chain]][2,,])))]
      if(any(is.nan(diag(x.ac[[chain]][2,,]))))
        cat(paste("WARNING: Variables",
                  paste(didnt.mix,collapse=", "),
                  "did not mix AT ALL. MCMC should be rerun with different proposal parameters!\n"))
      for(i in 1:2){
        cat(paste(dimnames(x.ac[[chain]])[[1]][i],"\n"))
        print(x.ac[[chain]][i,,])
        cat("\n")
      }
    }
  }

  if("acf" %in% which.diags)
    autocorr.plot(x)

  if("trace" %in% which.diags){
    plot(x)
  }

  if("raftery" %in% which.diags){
    #' @importFrom coda raftery.diag
    rd<-try(raftery.diag(x,r=0.0125))
    if(inherits(rd,"try-error")){
      cat("Raftery-Lewis diagnostic failed, likely due to some of the vairables not mixing at all.\n MCMC should be rerun.\n")
      return(invisible(NULL))
    }
    print(rd)
    invisible(rd)
  }
}

#' Convert an ERGMM Object to an MCMC list object for Diagnostics.
#' 
#' Functions to extract a subset of MCMC-sampled variables from an object of
#' class \code{\link[=ergmm.object]{ergmm}} and construct an
#' \code{\link[coda]{mcmc.list}} object.
#' 
#' Unless \code{which.vars} is specified, the \code{\link[coda]{mcmc.list}}
#' returned also includes all of the covariate coefficients.
#' 
#' Regardless of whether the MCMC run was single- or multi-threaded, this
#' function returns an \code{\link[coda]{mcmc.list}}, with a single thread, if
#' necessary.
#' 
#' @aliases as.mcmc.ergmm as.mcmc.list.ergmm
#' @param x An object of class \code{\link[=ergmm.object]{ergmm}}.
#' @param burnin If \code{TRUE}, generates an \code{\link[coda]{mcmc.list}}
#' object for the burnin (if stored) instead of the main sampling run.
#' @param which.vars A named list mapping variable names to the indices to
#' include. If given, overrides the defaults and all arguments that follow.
#' @param vertex.i A numeric vector of vertices whose latent space coordinates
#' and random effects to include.
#' @param ... Not used at this time.
#' @return A \code{\link[coda]{mcmc.list}} object with the sample of the
#' selected subset of the variables.
#' @seealso \code{\link{ergmm}}, \code{\link[coda]{mcmc.list}},
#' \code{\link{mcmc.diagnostics.ergmm}}
#' @keywords graphs debugging manip
#' @examples
#' 
#' \donttest{
#' library(coda)
#' data(sampson)
#' monks.fit<-ergmm(samplike~euclidean(d=2,G=3))
#' monks.fit.mcmc<-as.mcmc.list(monks.fit)
#' plot(monks.fit.mcmc)
#' raftery.diag(monks.fit.mcmc)
#' }
#' @importFrom coda mcmc as.mcmc mcmc.list as.mcmc.list
#' @importFrom ergm mcmc.diagnostics
#' @export
as.mcmc.list.ergmm<-function(x,burnin=FALSE,
                             which.vars=NULL,
                             vertex.i=c(1),...){
  extraneous.argcheck(...)
  n<-network.size(x[["model"]][["Yg"]])
  G<-x[["model"]][["G"]]
  d<-x[["model"]][["d"]]
  p<-x[["model"]][["p"]]
  start<-x[["control"]][["burnin"]]
  thin<-x[["control"]][["interval"]]

  as.mcmc.list.ergmm.par.list(if(burnin) x[["burnin.samples"]][[burnin]] else x[["sample"]],
                              if(is.null(which.vars)) list(lpY=1,
                                                           beta=seq_len(p),
                                                           Z=cbind(rep(vertex.i,each=d),rep(seq_len(d),length(vertex.i))),
                                                           sender=vertex.i,
                                                           receiver=vertex.i,
                                                           sociality=vertex.i,
                                                           dispersion=1) else which.vars,
                              start,thin)
}

#' @export
as.mcmc.ergmm <- as.mcmc.list.ergmm
