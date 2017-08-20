#  File R/summary.ergmm.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#' ERGMM Fit Summaries
#' 
#' \code{summary.ergmm} prodcues a summary of an
#' \code{\link[=ergmm.object]{ergmm}} object, including point estimates,
#' standard errors, and BIC calculation.
#'
#' @details
#' Note that BIC computed for the random effects models uses the same
#' formualtion as Handcock et al., so it is likely correct, but has not been
#' peer-reviewed.
#' 
#' This BIC can be (reasonably) safely used to select the number of clusters or
#' which fixed effects to include in the model. It is not clear whether it is
#' appropriate to use this BIC to select the dimension of latent space and
#' whether or not to include random actor effects. These considerations are
#' independent of the bug described below.
#' 
#' Prior to version 2.7.0, there was a bug in BIC calculation that used \eqn{p
#' + n(d+r+s)} as the number of parameters in the likelihood (where \eqn{p} is
#' the number of fixed effects, \eqn{n} the number of actors, \eqn{d}, the
#' latent space dimension, and \eqn{r} and \eqn{s} indicators of presence of
#' sender and receiver (or sociality) effects). This value should have been
#' just \eqn{p}.
#' 
#' The following applications could have produced different results:
#' \itemize{
#' \item{Using the BIC to select latent space dimension.}
#' \item{Using the BIC to decide whether or not to include random effects.}
#' }
#' The following applications could not (i.e., would be off by a constant):
#' \itemize{
#' \item{Using the BIC to select the number of clusters.}
#' \item{Using the BIC to select the fixed effects to be used.}
#' }
#' 
#' @aliases summary.ergmm print.summary.ergmm summary.ergmm.object bic.ergmm
#' @param object An \code{\link[=ergmm.object]{ergmm}} object to be summarized.
#' @param point.est Point estimates to compute: a character vector with some
#' subset of \code{"mle"}, \code{"pmean"}, \code{"mkl"}, and \code{"pmode"}.
#' Defaults to a concatenation of \code{"mle"} (if fit), \code{"pmean"}, and
#' \code{"mkl"} (if MCMC was run).
#' @param quantiles Posterior quantiles (credible intervals) to compute.
#' @param se Whether to compute standard errors. Defaults to \code{TRUE} if MLE
#' was fit.
#' @param eff.obs,bic.eff.obs What effective sample size to use for BIC
#' calculation?
#' \describe{
#'      \item{\code{"ties"}}{the number of non-missing ties in the network. This is the approach recommended by Handcock et al. (2007) and the default. Not well-defined for valued networks.}
#'      \item{\code{"dyads"}}{the number of non-missing dyads (potential ties) in the network.}
#'      \item{\code{"actors"}}{the number of actors in the network. The default prior to 2.7.0.}
#'      \item{a number}{to specify a specific sample size.}
#'      \item{\code{NULL}}{Don't compute the BIC at all. Mostly for internal use.}
#'    }
#' @param \dots Additional arguments.
#' @return For \code{summary}, an object of class
#' \code{\link[=summary.ergmm.object]{summary.ergmm}}. A print method is
#' available.
#' 
#' The BICs are available as the element "bic" of the object returned.
#' 
#' \code{bic.ergmm} returns the BIC for the model directly.
#' @seealso \code{\link{ergmm.object}}, \code{\link{ergmm}}
#' @references Chris Fraley and Adrian E. Raftery (2002). \emph{Model-based
#' clustering, discriminant analysis, and density estimation}. Journal of the
#' American Statistical Association, 97(458), 611-631.
#' 
#' Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum (2007).
#' \emph{Model-Based Clustering for Social Networks}.  Journal of the Royal
#' Statistical Society: Series A (Statistics in Society), 170(2), 301-354.
#' @keywords graphs models print
#' @examples
#' 
#' \donttest{
#' data(sampson)
#' # Fit the model for cluster sizes 1 through 4:
#' fits<-list(
#'            ergmm(samplike~euclidean(d=2,G=1)),
#'            ergmm(samplike~euclidean(d=2,G=2)),
#'            ergmm(samplike~euclidean(d=2,G=3)),
#'            ergmm(samplike~euclidean(d=2,G=4))
#'            )
#' 
#' \dontrun{
#' # Optionally, plot all fits.
#' lapply(fits,plot)
#' }
#' 
#' # Compute the BICs for the fits and plot them:
#' (bics<-reshape(
#'     as.data.frame(t(sapply(fits,
#'                            function(x)c(G=x$model$G,unlist(bic.ergmm(x))[c("Y","Z","overall")])))),
#'     list(c("Y","Z","overall")),idvar="G",v.names="BIC",timevar="Component",
#'     times=c("likelihood","clustering","overall"),direction="long"
#'     ))
#' 
#' with(bics,interaction.plot(G,Component,BIC,type="b",xlab="Clusters", ylab="BIC"))
#' 
#' # Summarize and plot whichever fit has the lowest overall BIC:
#' bestG<-with(bics[bics$Component=="overall",],G[which.min(BIC)])
#' summary(fits[[bestG]])
#' plot(fits[[bestG]])
#' }
#' 
#' @export
summary.ergmm <- function (object, point.est=c(
                                     if(!is.null(object[["mle"]])) "mle",
                                     if(!is.null(object[["sample"]])) c("pmean","mkl")),
                           quantiles=c(.025,.975),se="mle"%in%point.est,
                           bic.eff.obs=c("ties", "dyads", "actors"), ...)
{
  extraneous.argcheck(...)
  ## Just for convenience.
  sample<-object[["sample"]]
  model<-object[["model"]]
  control<-object[["control"]]
  d<-model[["d"]]
  p<-model[["p"]]
  G<-model[["G"]]
  n<-network.size(model[["Yg"]])
  sample.size<-control[["sample.size"]]
  

  summ<-list(ergmm=object,model=model)

  ## Compute the two-stage MLE point estimates.
  if("mle" %in% point.est){
    if(is.null(object[["mle"]])){
      stop("MLE was not computed for this fit.")
    }
    if(is.null(object[["mle"]][["cov"]])){
      if(model[["sender"]] || model[["receiver"]] || model[["sociality"]])
        warning("Fitting random effects as fixed effects.")
      if(se){
        ## Refit the MLE (mostly for the Hessian)
        mle<-find.mle(model,object[["mle"]],control=control,hessian=TRUE)
      }
      else mle<-object[["mle"]]

      if(d>0) mle[["Z"]]<-scale(mle[["Z"]],scale=FALSE)
      if(p>0){
        if(se){
          beta.hess<-mle[["hessian"]][1:p,1:p]
          
          beta.cov <- try(MASS::ginv(-beta.hess), silent=TRUE)
          if(inherits(beta.cov,"try-error")){
            warning("Coefficient Hessian appears to be singular. Using a less accurate estimate.")
            beta.cov <- diag(1/diag(-beta.hess))
          }
          
          colnames(beta.cov) <- rownames(beta.cov) <- model[["coef.names"]]
        
          mle[["cov"]]<-beta.cov
          #' @importFrom stats cov2cor
          mle[["cor"]]<-cov2cor(beta.cov)
        
          z.values<-mle[["beta"]]/sqrt(diag(mle[["cov"]]))
          #' @importFrom stats pnorm
          coef.table<-data.frame(mle[["beta"]],sqrt(diag(mle[["cov"]])),
                                 z.values,2*pnorm(abs(z.values),0,1,lower.tail=FALSE),
                                 row.names=model[["coef.names"]])
        }
        else coef.table<-data.frame(mle[["beta"]],
                                    row.names=model[["coef.names"]])
        colnames(coef.table)<-c("Estimate",if(se) "Std. Error",if(se) "z value",if(se) "Pr(>|z|)")
        mle[["coef.table"]]<-coef.table
      }

      if(model[["sender"]]){
        if(model[["intercept"]]){
          mle[["beta"]][1]<-mle[["beta"]][1]-mean(mle[["sender"]])
          mle[["sender"]]<-mle[["sender"]]-mean(mle[["sender"]])
        }
        mle[["sender.var"]]<-mean(mle[["sender"]]^2)
      }
        
      if(model[["receiver"]]){
        if(model[["intercept"]]){
          mle[["beta"]][1]<-mle[["beta"]][1]-mean(mle[["receiver"]])
          mle[["receiver"]]<-mle[["receiver"]]-mean(mle[["receiver"]])
        }
        mle[["receiver.var"]]<-mean(mle[["receiver"]]^2)
      }

      if(model[["sociality"]]){
        if(model[["intercept"]]){
          mle[["beta"]][1]<-mle[["beta"]][1]-mean(mle[["sociality"]])
          mle[["sociality"]]<-mle[["sociality"]]-mean(mle[["sociality"]])
        }
        mle[["sociality.var"]]<-mean(mle[["sociality"]]^2)
      }

      if(model[["dispersion"]]){
        mle[["dispersion"]]<-mle[["dispersion"]]
      }
      
      object[["mle"]]<-mle
    }
    summ[["mle"]]<-object[["mle"]]
  }

  ## Compute the posterior mean point estimates.
  if("pmean" %in% point.est){
    if(is.null(object[["sample"]])){
      stop("MCMC was not was not run for this fit.")
    }
    if(is.null(object[["pmean"]])){
      pmean<-list()
      for(name in names(sample)){
        if(is.null(sample[[name]])) next
        name.dim<-length(dim(sample[[name]]))
        pmean[[name]]<-{
          if(name.dim<2) mean(sample[[name]])
          else if(name.dim==2) apply(sample[[name]],2,mean)
          else if(name.dim==3) apply(sample[[name]],2:3,mean)
        }
      }
      if(G>0){
        pmean[["Z.pZK"]]<-t(apply(sample[["Z.K"]],2,tabulate,G))/sample.size
        pmean[["Z.K"]]<-apply(pmean[["Z.pZK"]],1,which.max)
      }
    
      beta.cov<-cov(sample[["beta"]])
      colnames(beta.cov) <- rownames(beta.cov) <- model[["coef.names"]]
      
      pmean[["cov"]]<-beta.cov
      pmean[["cor"]]<-cov2cor(beta.cov)
      
      beta.q0<-apply(sample[["beta"]],2,function(x) min(mean(x<=0),mean(x>=0))*2)
      
      coef.table<-data.frame(pmean[["beta"]],
                           t(apply(sample[["beta"]],2,function(x)quantile(x,quantiles))),
                           beta.q0,
                           row.names=model[["coef.names"]])
      colnames(coef.table)<-c("Estimate",paste(quantiles*100,"%",sep=""),"2*min(Pr(>0),Pr(<0))")
      pmean[["coef.table"]]<-coef.table
      object[["pmean"]]<-pmean
    }
    summ[["pmean"]]<-object[["pmean"]]
  }
  ## Compute the MKL point estimates.
  if("mkl" %in% point.est){
    if(is.null(object[["mkl"]])){
      stop("MKL was not produced for this fit.")
    }
    mkl<-summ[["mkl"]]<-object[["mkl"]]
    coef.table<-data.frame(mkl[["beta"]],
                           row.names=model[["coef.names"]])
    colnames(coef.table)<-c("Estimate")
    summ[["mkl"]][["coef.table"]]<-coef.table
  }

  if("pmode" %in% point.est){
    if(is.null(object[["pmode"]])){
      stop("Conditional posterior mode was not computed for this fit.")
    }
    summ[["pmode"]]<-object[["pmode"]]
  }

  if(!is.null(object[["mkl"]]) && !is.null(bic.eff.obs)){
    summ[["bic"]]<-bic.ergmm(object, eff.obs = bic.eff.obs, ...)
  }

  class(summ)<-'summary.ergmm'
  summ
}

#' @export
print.summary.ergmm<-function(x,...){
  ## For convenience
  model<-x[["model"]]
  control<-x[["ergmm"]][["control"]]
  
  cat("\n==========================\n")
  cat("Summary of model fit\n")
  cat("==========================\n\n")

  cat("Formula:   ")
  print(model[["formula"]])
  cat("Attribute: ")
  if(is.null(model[["response"]])) cat("edges") else cat(model[["response"]])
  cat("\n")
  cat("Model:    ",model[["family"]],"\n")
  
  digits = max(3, getOption("digits") - 3)
  
  if(!is.null(x[["pmean"]])) cat ("MCMC sample of size ", control[["sample.size"]], ", draws are ",
       control[["interval"]]," iterations apart, after burnin of ",control[["burnin"]], " iterations.\n",sep="")
       
  if(!is.null(x[["pmean"]])){
    cat("Covariate coefficients posterior means:\n")
    #' @importFrom stats printCoefmat
    printCoefmat(as.matrix(x[["pmean"]][["coef.table"]]),P.values=TRUE,has.Pvalue=TRUE)
    cat("\n")
    if(!is.null(x[["pmean"]][["dispersion"]]))
      cat("Dispersion parameter: ",x[["pmean"]][["dispersion"]],".\n", sep="")
    if(!is.null(x[["pmean"]][["sender.var"]]))
      cat("Sender effect variance: ",x[["pmean"]][["sender.var"]],".\n", sep="")
    if(!is.null(x[["pmean"]][["receiver.var"]]))
      cat("Receiver effect variance: ",x[["pmean"]][["receiver.var"]],".\n", sep="")
    if(!is.null(x[["pmean"]][["sociality.var"]]))
      cat("Sociality effect variance: ",x[["pmean"]][["sociality.var"]],".\n", sep="")
  }

  if(!is.null(x[["mle"]])){
    cat("Covariate coefficients MLE:\n")
    printCoefmat(as.matrix(x[["mle"]][["coef.table"]]),P.values=length(names(x[["mle"]][["coef.table"]]))>1)
    cat("\n")
    if(!is.null(x[["mle"]][["dispersion"]]))
      cat("Dispersion parameter: ",x[["mle"]][["dispersion"]],".\n", sep="")
  }
  if(!is.null(x[["bic"]])){
    cat("Overall BIC:       ", x[["bic"]][["overall"]],"\n")
    cat("Likelihood BIC:    ", x[["bic"]][["Y"]],"\n")
    if(model[["d"]]>0){
      cat("Latent space/clustering BIC:    ", x[["bic"]][["Z"]],"\n")
    }
    if(model[["sender"]]){
      cat("Sender effect BIC:    ", x[["bic"]][["sender"]],"\n")
    }
    if(model[["receiver"]]){
      cat("Receiver effect BIC:    ", x[["bic"]][["receiver"]],"\n")
    }
    if(model[["sociality"]]){
      cat("Sociality effect BIC:    ", x[["bic"]][["sociality"]],"\n")
    }
    cat("\n")
  }

  if(!is.null(x[["mkl"]])){
    cat("Covariate coefficients MKL:\n")
    print(x[["mkl"]][["coef.table"]])
    cat("\n")
    if(!is.null(x[["mkl"]][["dispersion"]]))
      cat("Dispersion parameter: ",x[["mkl"]][["dispersion"]]," (probably invalid).\n", sep="")
    cat("\n")
  }
  if(!is.null(x[["pmode"]])){
    cat("Covariate coefficients posterior mode:\n")
    print(x[["pmode"]][["coef.table"]])
    cat("\n")
    if(!is.null(x[["pmode"]][["dispersion"]]))
      cat("Dispersion parameter: ",x[["pmode"]][["dispersion"]],".\n", sep="")
    cat("\n")
  }
}

#' @export bic.ergmm
bic.ergmm<-function(object, eff.obs=c("ties", "dyads", "actors"), ...){
  extraneous.argcheck(...)
  if(is.null(object[["mkl"]])){
    stop("MKL estimates were not computed for this fit.")
  }

  model<-object[["model"]]

  n<-network.size(model[["Yg"]])
  
  if(is.character(eff.obs)){
    eff.obs <- switch(match.arg(eff.obs),
                      ties = {if(!all(match(model[["Ym"]], c(0,1)) > 0, na.rm=TRUE)) warning('Number of "ties" in a valued network may not be well-defined.'); network.edgecount(model[["Yg"]])},
                      dyads = network.dyadcount(model[["Yg"]]),
                      actors = n)
  }
  
  
  condZRE<-with(object,find.mle(model,mkl,given=list(Z=mkl[["Z"]],sender=mkl[["sender"]],receiver=mkl[["receiver"]],sociality=mkl[["sociality"]]),control=object[["control"]]))

  bic<-with(model,list(Y = -2*condZRE[["lpY"]] + (p)*log(eff.obs),
                              Z =
                              if(d>0){
                                if(G>0){
                                  mbc.llk<-Inf
                                  Gsub<--1
                                  while(!is.finite(mbc.llk)){
                                    Gsub<-Gsub+1
                                    mbc.llk<-mbc.VII.EM(G-Gsub,object[["mkl"]][["Z"]])[["llk"]]
                                  }
                                  if(Gsub) warning(paste("Bad clustering: treating",Gsub,"clusters as empty."))
                                  -2*mbc.llk+((G-Gsub)-1 + d*(G-Gsub) + (G-Gsub))*log(n)
                                } else {
                                  -2*sum(dnorm(object[["mkl"]][["Z"]],0,sqrt(mean(condZRE[["Z"]]^2)*d),log=TRUE))+1*log(n*d)
                                }
                              } else 0,
                              sender=if(sender) -2*sum(dnorm(condZRE[["sender"]],0,sqrt(mean(condZRE[["sender"]]^2)),log=TRUE))+1*log(n) else 0,
                              receiver=if(receiver) -2*sum(dnorm(condZRE[["receiver"]],0,sqrt(mean(condZRE[["receiver"]]^2)),log=TRUE))+1*log(n) else 0,
                              sociality=if(sociality) -2*sum(dnorm(condZRE[["sociality"]],0,sqrt(mean(condZRE[["sociality"]]^2)),log=TRUE))+1*log(n) else 0
                              )
            )

  if(!.latentnetEnv$BIC.warned){
    message("NOTE: It is not certain whether it is appropriate to use latentnet's BIC to select latent space dimension, whether or not to include actor-specific random effects, and to compare clustered models with the unclustered model.")
    .latentnetEnv$BIC.warned <- TRUE
  }
  
  bic[["overall"]]<-sum(unlist(bic))
  
  bic
}
