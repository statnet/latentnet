#  File R/ergmm.par.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
# Utilities for dealing with MCMC output produced by *.MCMC.C functions.

ERGMM.PAR_VAR_NAMES<-c("beta","Z","sender","receiver","sociality",
                       "Z.var","Z.mean","Z.K",
                       "sender.var","receiver.var","sociality.var",
                       "dispersion")
ERGMM.PAR_LLK_NAMES<-c("beta","Z","sender","receiver","sociality","dispersion")

del.iteration<-function(mcmcsample,i){
  for(name in names(mcmcsample)){
    if(length(mcmcsample[[name]])>0){
      if(length(dim(mcmcsample[[name]]))<=1) mcmcsample[[name]]<-mcmcsample[[name]][-i]
      else if(length(dim(mcmcsample[[name]]))==2) mcmcsample[[name]]<-mcmcsample[[name]][-i,,drop=FALSE]
      else if(length(dim(mcmcsample[[name]]))==3) mcmcsample[[name]]<-mcmcsample[[name]][-i,,,drop=FALSE]
    }
  }
  mcmcsample
}

seldrop<-function(x,i){
  array(c(x),dim=dim(x)[-i])
}

#' @name ergmm.par.list
#'
#' @title A List of ERGMM Parameter Configuration
#'
#' @description
#' A class \code{\link[=ergmm.par.list.object]{ergmm.par.list}} to represent a
#' series of parameter configurations for the same exponential random graph
#' mixed model.
#'
#' @param x an \code{ergm.par.list} object.
#' @param i index of the iteration to extract.
#' @param ... extra arguments, currently unused.
#' 
#' @details
#'  \code{[[} operator with a
#'  numeric or integer index returns a
#'  list with the the
#'  configuration with that index. \code{[} operator given a numeric
#'  vector returns a \code{ergmm.par.list} object with the subset of
#'  configurations with the indices given.
#'  
#'  The structure of \code{ergmm.par.list} is derived from named lists, with each entry having an
#'  additional dimension (always the first one), indexed by
#'  configuration. That is, scalars become vectors, vectors become
#'  matrixes with the original vectors in rows, and matrixes become
#'  3-dimensional arrays, with the original matrixes indexed by their
#'  first dimension. See \code{\link{terms.ergmm}} for comon elements of
#'  these configurations.
#'  In some cases, such as when representing MCMC or optimization output,
#'  the object may also have some of the following elements:
#'  \describe{
#'    \item{\code{mlp}}{\eqn{\log p(Y,Z,\beta,\mu,\sigma,\delta,\gamma,\sigma_\delta,\sigma_\gamma,|K)},
#'      the joint
#'      probability/density of network, the covariate coefficients, the
#'      latent space positions and parameters, and the random effects and
#'      their variances, conditional on cluster assignments.}
#'    \item{\code{lpY}}{\eqn{\log p(Y|\dots)}, depending on the model, the log-probability or
#'      log-density of the network conditional on all the parameters.}
#'    \item{\code{lpZ}}{\eqn{\log p(Z|\mu,\sigma,K)}, the log-density of latent space positions conditional on
#'      latent space or cluster parameters and cluster assignments.}
#'    \item{\code{lpbeta}}{\eqn{\log p(\beta)}, the prior log-density of
#'      the covariate coefficients.}
#'    \item{\code{lpRE}}{\eqn{\log p(\delta,\gamma|\sigma_\delta,\sigma_\gamma)}, the log-density of all random effects, conditional on
#'      their respective variances.}
#'    \item{\code{lpLV}}{\eqn{\log p(\mu,\sigma)}, the prior log-density
#'      of latent space or cluster parameters (but not that of the cluster
#'      assignments).}
#'    \item{\code{lpREV}}{\eqn{\log p(\sigma_\delta,\sigma_\gamma)}, the prior log-density of all random effect variances.}
#'    \item{\code{Z.rate}}{Proportion of single-vertex proposals accepted over the preceding
#'      interval.}
#'    \item{\code{beta.rate}}{
#'      Proportion of group proposals accepted over the preceding interval.}
#'  }
#' @aliases ergmm.par.list.object ergmm.par.list as.ergmm.par.list
#' as.mcmc.list.ergmm.par.list [.ergmm.par.list [[.ergmm.par.list
#' $.ergmm.par.list length.ergmm.par.list unstack.ergmm.par.list del.iteration
#' @seealso \code{\link{ergmm}}
#' @keywords graphs utilities methods list manip
NULL

#' @rdname ergmm.par.list
#' @export
as.ergmm.par.list<-function(x,...){
  class(x)<-"ergmm.par.list"
  x
}

#' @rdname ergmm.par.list
#' @export
"[.ergmm.par.list"<-function(x,i){
  if(!is.numeric(i)) stop("Index vector to the '[' operator of an ergmm.par.list must be of mode numeric or integer.")
  l<-list()
  
  for(name in names(x)){
    if(length(x[[name]])>0){
      d<-dim(x[[name]])
      if(length(d)<=1) l[[name]]<-x[[name]][i]
      else if(length(d)==2) l[[name]]<-x[[name]][i,,drop=FALSE]
      else if(length(d)==3) l[[name]]<-x[[name]][i,,,drop=FALSE]
    }
  }
  class(l)<-"ergmm.par.list"
  l
}

#' @rdname ergmm.par.list
#' @export
length.ergmm.par.list<-function(x){
  if(is.null(dim(x[[names(x)[1]]]))) length(x[[names(x)[1]]])
  else dim(x[[names(x)[1]]])[1]
}

#' @rdname ergmm.par.list
#' @export
`[[.ergmm.par.list`<-`$.ergmm.par.list`<-function(x,i){
  ## Delete its class, to keep it from recursing.
  tmp<-class(x)
  class(x)<-NULL
  if(inherits(i,"character")){ ## If the index is a character, return all the draws for the corresponding variable.
    xi<-x[i][[1]]
    class(x)<-tmp
    return(xi)
  }
  else{  ## If it's a number, return a configuration with that iteration number.
    l<-list()
    ## Do NOT seldrop 1D parameters. Those actually need to become vectors.
    ## (As opposed to becoming 1*n matrices.)
    
    for(name in names(x)){
      if(length(x[[name]])>0){
        if(length(dim(x[[name]]))<=1) l[[name]]<-x[[name]][i]
        else if(length(dim(x[[name]]))==2) l[[name]]<-x[[name]][i,]
        else if(length(dim(x[[name]]))==3) l[[name]]<-seldrop(x[[name]][i,,,drop=FALSE],1)
      }
    }
    class(x)<-tmp
    return(l)
  }
}

.stack.ergmm.par.list.list<-function(x,...){
  extraneous.argcheck(...)
  mcmcsample<-list()

  for(name in names(x[[1]]))
    mcmcsample[[name]]<-abind::abind(sapply(seq_along(x),
                                            function(i) x[[i]][[name]],
                                            simplify=FALSE),along=1)

  attr(mcmcsample,"breaks")<-cumsum(c(sapply(seq_along(x),
                                             function(i) length(x[[i]]),
                                             simplify=FALSE)))
  class(mcmcsample)<-"ergmm.par.list"
  mcmcsample
}

#' @importFrom utils unstack
#' @rdname ergmm.par.list
#' @export
unstack.ergmm.par.list<-function(x,...){
  extraneous.argcheck(...)
  mcmcList<-list()

  if(is.null(attr(x,"breaks"))){
    mcmcList[[1]]<-x
  }
  else{  
    breaks<-c(0,attr(x,"breaks"))
    
    for(i in seq_along(breaks[-1])){
      mcmcList[[i]]<-x[(breaks[i]+1):breaks[i+1]]
    }
  }
  mcmcList
}

as.mcmc.list.ergmm.par.list<-function(x,which.vars,start=1,thin=1,...){
  x<-unstack(x)
  m.l<-list()
  for(thread in seq_along(x)){
    S<-length(x[[thread]])
    m<-matrix(numeric(0),S,0)
    for(name in names(which.vars)){
      if(length(x[[thread]][[name]])==0) next
      if(length(dim(x[[thread]][[name]]))<=1){
        m2<-cbind(x[[thread]][[name]])
        colnames(m2)<-name
        m<-cbind(m,m2)
      }else if(length(dim(x[[thread]][[name]]))==2){
        m2<-x[[thread]][[name]][,c(which.vars[[name]]),drop=FALSE]
        colnames(m2)<-paste(name,vapply(which.vars[[name]],function(x) paste(x,sep='.'), character(1)),sep='.')
        m<-cbind(m,m2)
      }else if(length(dim(x[[thread]][[name]]))==3){
        for(i in 1:dim(which.vars[[name]])[1]){
          i<-which.vars[[name]][i,]
          m2<-cbind(x[[thread]][[name]][,i[1],i[2]])
          colnames(m2)<-paste(name,paste(i,sep='.',collapse='.'),sep='.')
          m<-cbind(m,m2)
        }
      }
    }
    m<-mcmc(m,start=start,thin=thin)
    m.l[[thread]]<-m
  }
  eval(as.call(c(mcmc.list,m.l)))
}
