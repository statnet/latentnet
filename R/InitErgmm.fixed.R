#  File R/InitErgmm.fixed.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

#' @importFrom statnet.common ult
nonlatent_error <- function(...){
  myname <- sub("InitErgmTerm.", "", ult(as.character(sys.call()[[1]])), fixed=TRUE)

  stop("It appears that you are using ", sQuote("latentnet"), " term ", sQuote(myname)," in a non-latent-space context. They should only be used in ", sQuote("ergmm()"), " calls and similar.", call. = FALSE)
}

.ergmm.add.fixed<-function(model, X, mean, var, coef.names=NULL, where=c("append","prepend")){
  where <- match.arg(where)
  
  if(length(dim(X))==2) dim(X) <- c(dim(X),1)
  if(!is.null(coef.names)) dimnames(X) <- list(NULL, NULL, coef.names)

  p <- dim(X)[3]
  model[["p"]]<-model[["p"]]+p
  model[["coef.names"]]<-switch(where,
                                append = c(model[["coef.names"]], dimnames(X)[[3]]),
                                prepend = c(dimnames(X)[[3]], model[["coef.names"]])
                                )
                                
  model[["prior"]][["beta.mean"]]<-switch(where,
                                          append = c(model[["prior"]][["beta.mean"]], rep(mean,length.out=p)),
                                          prepend = c(rep(mean,length.out=p), model[["prior"]][["beta.mean"]])
                                          )
  
  model[["prior"]][["beta.var"]]<-switch(where,
                                         append = c(model[["prior"]][["beta.var"]], rep(var,length.out=p)),
                                         prepend = c(rep(var,length.out=p), model[["prior"]][["beta.var"]])
                                         )
                                         

  Yg <- model[["Yg"]]


  Xtmp <- list()
  for(i in seq_len(p)){
    xm <- seldrop(X[,,i,drop=FALSE],3)
    # If the network is undirected, symmetrize:
    if(!is.directed(Yg)) xm[lower.tri(xm)] <- t(xm)[lower.tri(xm)]
    
    # If the network is bipartite and the matrix has b1*b2 dimensions,
    # it needs to be augmented to (b1+b2)*(b1+b2):
    if(is.bipartite(Yg)
       && all(dim(xm)==c(Yg%n%"bipartite", network.size(Yg)-Yg%n%"bipartite")))
      xm <- bipartite.augment(xm)
    
    Xtmp <- c(Xtmp,list(xm))
  }

  model[["X"]]<-switch(where,
                       append = c(model[["X"]],Xtmp),
                       prepend = c(Xtmp, model[["X"]])
                       )

  model
}

.import.ergm.term <-function(model, term.index, term.name, ..., mean=0, var=9){
  Yg<-model[["Yg"]]
  f <- ~Yg

  #' @importFrom statnet.common list_rhs.formula
  f[[3]] <- list_rhs.formula(model$formula)[[term.index]]

  if(!is.dyad.independent(f)) warning("Term `", term.name, "` induces dyadic dependence. Likelihood will be effectively replaced by pseudolikelihood.", call.=FALSE)
  if(has.loops(Yg)) warning("Imported ergm term `", term.name, "` will set its dyadic covariate for self-loops, X[i,i,k], to 0. Use `loopfactor` and `loopcov` to model self-loops.", call.=FALSE)
  
  X <- ergmMPLE(f, output="array")$predictor
  
  X[c(is.na(X))] <- 0

  .ergmm.add.fixed(model, X, mean, var)
}

#' @templateVar name Intercept
#' @aliases 1-ergmTerm intercept-ergmTerm
#' @title Intercept
#' @description This term serves as an intercept term, is included by
#'   default (though, as in \code{\link{lm}}, it can be excluded by
#'   adding \code{+0} or \code{-1} into the model formula). It adds
#'   one covariate to the model, for which \code{x[i,j]=1} for all
#'   \code{i} and \code{j}.
#'
#'   It can be used explicitly to set prior mean and variance for the
#'   intercept term.
#'
#'   This term differs from the \code{ergm}'s
#'   \code{\link{edges-ergmTerm}} term if the network has self-loops.
#'
#' @usage
#' # binary: 1(mean=0, var=9)
#'
#' # binary: Intercept(mean=0, var=9)
#'
#' # binary: intercept(mean=0, var=9)
#'
#' # valued: 1(mean=0, var=9)
#'
#' # valued: Intercept(mean=0, var=9)
#'
#' # valued: intercept(mean=0, var=9)
#'
#' @template ergmTerm-latentnet-prior
#'
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @concept latent
InitErgmTerm.Intercept <- InitErgmTerm.intercept <- InitErgmTerm.1 <- nonlatent_error
InitErgmm.Intercept<-InitErgmm.intercept<-InitErgmm.1<-function(model, mean=0, var=9){
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 1:4))
    stop(paste("`edges` model term expected between 0 and 2 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  .ergmm.add.fixed(model,
                   matrix(1,network.size(model[["Yg"]]),network.size(model[["Yg"]])),
                   mean, var,
                   "(Intercept)")
}

#' @templateVar name loops
#' @title Self-loops
#' @description Effect of the dyad being a self-loop (i.e., \eqn{(i,i)}).
#'
#' @usage
#' # binary: loops(mean=0, var=9)
#'
#' # valued: loops(mean=0, var=9)
#'
#' @template ergmTerm-latentnet-prior
#'
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @concept latent
InitErgmTerm.loops <- nonlatent_error
InitErgmm.loops<-function (model, mean=0, var=9){
  if(!has.loops(model[["Yg"]]))
    stop("Self-loop term is  meaningless in a network without self-loops", call.=FALSE)

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 1:3))
    stop(paste("`loops` model term expected between 0 and 2 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  .ergmm.add.fixed(model, diag(1,network.size(model[["Yg"]]),network.size(model[["Yg"]])), mean, var, "loops")
}

#' @templateVar name loopcov
#' @title Covariate effect on self-loops
#' @description This term adds one covariate to the model, for which
#'   \code{x[i,i]=attrname(i)} and \code{x[i,j]=0} for \code{i!=j}.
#'   This term only makes sense if the network has self-loops.
#'
#' @usage
#' # binary: loopcov(attrname, mean=0, var=9)
#'
#' # valued: loopcov(attrname, mean=0, var=9)
#'
#' @param attrname a character string giving the name of a numeric
#'   (not categorical) attribute in the network's vertex attribute
#'   list.
#'
#' @template ergmTerm-latentnet-prior
#'
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#'
#' @concept dyad-independent
#' @concept quantitative nodal attribute
#' @concept directed
#' @concept undirected
#' @concept latent
InitErgmTerm.loopcov <- nonlatent_error
InitErgmm.loopcov <- function (model, attrname, mean=0, var=9){
  if(!has.loops(model[["Yg"]]))
    stop("Self-loop covariates are meaningless in a network without self-loops", call.=FALSE)

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:4))
    stop(paste("loopcov() model term expected between 1 and 3 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname

  xm<-diag(x,n,n)
  cn<-paste("loopcov",attrname,sep=".")

  .ergmm.add.fixed(model, xm, mean, var, cn)
}

#' @templateVar name loopfactor
#' @title Factor attribute effect on self-loops
#' @description This term adds multiple covariates to the model, one
#'   for each of (a subset of) the unique values of the
#'   \code{attrname} attribute (or each combination of the attributes
#'   given). Each of these covariates has \code{x[i,i]=1} if
#'   \code{attrname(i)==l}, where \code{l} is that covariate's level,
#'   and \code{x[i,j]=0} otherwise. To include all attribute values se
#'   \code{base=0} -- because the sum of all such statistics equals
#'   twice the number of self-loops and hence a linear dependency
#'   would arise in any model also including \code{loops}. Thus, the
#'   \code{base} argument tells which value(s) (numbered in order
#'   according to the \code{sort} function) should be omitted. The
#'   default value, \code{base=1}, means that the smallest (i.e.,
#'   first in sorted order) attribute value is omitted. For example,
#'   if the \dQuote{fruit} factor has levels \dQuote{orange},
#'   \dQuote{apple}, \dQuote{banana}, and \dQuote{pear}, then to add
#'   just two terms, one for \dQuote{apple} and one for \dQuote{pear},
#'   then set \dQuote{banana} and \dQuote{orange} to the base
#'   (remember to sort the values first) by using
#'   \code{nodefactor("fruit", base=2:3)}. For an analogous term for
#'   quantitative vertex attributes, see
#'   \code{nodecov}.\code{attrname} is a character string giving the
#'   name of a numeric (not categorical) attribute in the network's
#'   vertex attribute list.  This term adds one covariate to the
#'   model, for which \code{x[i,i]=attrname(i)} and \code{x[i,j]=0}
#'   for \code{i!=j}.  This term only makes sense if the network has
#'   self-loops.
#'
#' @usage
#' # binary: loopfactor(attrname, mean=0, var=9)
#'
#' # valued: loopfactor(attrname, mean=0, var=9)
#'
#' @param attrname a character vector giving one or more names of
#'   categorical attributes in the network's vertex attribute list.
#'
#' @template ergmTerm-latentnet-prior
#'
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#'
#' @concept dyad-independent
#' @concept categorical nodal attribute
#' @concept directed
#' @concept undirected
#' @concept latent
InitErgmTerm.loopfactor <- nonlatent_error
InitErgmm.loopfactor <- function (model, attrname, base=1, mean=0, var=9){
  if(!has.loops(model[["Yg"]]))
    stop("Self-loop covariates are meaningless in a network without self-loops", call.=FALSE)

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:4))
    stop(paste("loopfactor() model term expected between 1 and 3 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  n<-network.size(model[["Yg"]])

  if(length(attrname)==1)
    x<-model[["Yg"]] %v% attrname
  else
    do.call(paste, c(lapply(attrname, function(a) get.vertex.attribute(model[["Yg"]], a)), sep="."))

  ls<-sort(unique(x))
  if(NVL(base,0)!=0){
    ls <- ls[-base]
    if(length(ls)==0){
      warning("`loopfactor` term deleted because contributes no statistics.")
      return(model)
    }
  }
  
  mean<-rep(mean,length.out=length(ls))
  var<-rep(var,length.out=length(ls))
  for(li in seq_along(ls)){
    l<-ls[[li]]
    xm<-diag(x==l,n,n)
    cn<-paste('loopfactor',paste(attrname,collapse="."),l,sep=".")
    model <- .ergmm.add.fixed(model, xm, mean[li], var[li], cn)
  }
  model
}

#' @templateVar name latentcov
#' @title Edge covariates for the latent model
#' @description This term adds multiple covariates to the model, one
#'   for each of (a subset of) the unique values of the
#'   \code{attrname} attribute (or each combination of the attributes
#'   given). Each of these covariates has \code{x[i,i]=1} if
#'   \code{attrname(i)==l}, where \code{l} is that covariate's level,
#'   and \code{x[i,j]=0} otherwise. To include all attribute values se
#'   \code{base=0} -- because the sum of all such statistics equals
#'   twice the number of self-loops and hence a linear dependency
#'   would arise in any model also including \code{loops}. Thus, the
#'   \code{base} argument tells which value(s) (numbered in order
#'   according to the \code{sort} function) should be omitted. The
#'   default value, \code{base=1}, means that the smallest (i.e.,
#'   first in sorted order) attribute value is omitted. For example,
#'   if the \dQuote{fruit} factor has levels \dQuote{orange},
#'   \dQuote{apple}, \dQuote{banana}, and \dQuote{pear}, then to add
#'   just two terms, one for \dQuote{apple} and one for \dQuote{pear},
#'   then set \dQuote{banana} and \dQuote{orange} to the base
#'   (remember to sort the values first) by using
#'   \code{nodefactor("fruit", base=2:3)}. For an analogous term for
#'   quantitative vertex attributes, see
#'   \code{nodecov}.\code{attrname} is a character string giving the
#'   name of a numeric (not categorical) attribute in the network's
#'   vertex attribute list.  This term adds one covariate to the
#'   model, for which \code{x[i,i]=attrname(i)} and \code{x[i,j]=0}
#'   for \code{i!=j}.  This term only makes sense if the network has
#'   self-loops.
#'
#'   \code{latentcov} can be called more than once, to model the
#'   effects of multiple covariates. Note that some covariates can be
#'   more conveniently specified using the following terms.
#'
#' @usage
#' # binary: latentcov(x, attrname=NULL, mean=0, var=9)
#'
#' # valued: latentcov(x, attrname=NULL, mean=0, var=9)
#'
#' @param x either a matrix of covariates on each pair of vertices, a
#'   network, or an edge attribute.
#' @param attrname optional argument to provide the name of the edge attribute.
#'
#' @template ergmTerm-latentnet-prior
#'
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @concept latent
InitErgmTerm.latentcov <- nonlatent_error
InitErgmm.latentcov<-function (model, x, attrname=NULL,
                               mean=0, var=9) 
{
  if(!has.loops(model[["Yg"]])) warning("Term `latentcov` is deprecated for networks without self-loops. Use `edgecov` from package `ergm` instead.")
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("latentcov() model term expected between 1 and 4 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  #Coerce x to an adjacency matrix
  if(is.network(x)){
    xm<-as.matrix.network(x,matrix.type="adjacency",attrname)
    cn<-if(!is.null(attrname)) attrname else paste("network",length(model[["X"]])+1)
  }else if(is.character(x)){
    xm<-as.matrix.network(model[["Yg"]],matrix.type="adjacency",x)
    cn<-x
  }else{
    xm<-as.matrix(x)
    cn<-if(!is.null(attrname)) attrname else 
				paste("latentcov", as.character(sys.call(0)[[3]]),	sep = ".")
  }

  .ergmm.add.fixed(model, xm, mean, var, cn)
}

#' @templateVar name sendercov
#' @title Sender covariate effect
#' @description \emph{Deprecated for networks without self-loops. Use
#'   \code{\link{nodeocov-ergmTerm}},
#'   \code{\link{nodeofactor-ergmTerm}},
#'   \code{\link{nodecov-ergmTerm}} or
#'   \code{\link{nodefactor-ergmTerm}} instead.}
#'
#'   If the attribute is numeric, this term adds one covariate to the
#'   model equaling \code{attrname(i)}. If the attribute is not
#'   numeric or \code{force.factor==TRUE}, this term adds \eqn{p-1}
#'   covariates to the model, where \eqn{p} is the number of unique
#'   values of \code{attrname}.  The \eqn{k}th such covariate has the
#'   value \code{attrname(i) == value(k+1)}, where \code{value(k)} is
#'   the \eqn{k}th smallest unique value of the \code{attrname}
#'   attribute. This term only makes sense if the network is directed.
#'
#' @usage
#' # binary: sendercov(attrname, force.factor=FALSE, mean=0, var=9)
#'
#' # valued: sendercov(attrname, force.factor=FALSE, mean=0, var=9)
#'
#' @param attrname a character string giving the name of an attribute
#'   in the network's vertex attribute list.
#' @param force.factor logical, indicating if `attrname`'s value
#'   should be interpreted as categorical even if numeric.
#'
#' @template ergmTerm-latentnet-prior
#'
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#' @template ergmTerm-directed
#'
#' @concept dyad-independent
#' @concept directed
#' @concept latent
#' @concept loops
InitErgmTerm.sendercov <- nonlatent_error
InitErgmm.sendercov<-function (model, attrname, force.factor=FALSE, mean=0, var=9) 
{
  if(!has.loops(model[["Yg"]])) warning("Term `sendercov` is deprecated for networks without self-loops. Use `nodeocov`, `nodecov`, `nodeofactor`, or `nodefactor` from package `ergm` instead.")

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("sendercov() model term expected between 1 and 4 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  if (!is.directed(model[["Yg"]]))
    stop("Sender covariates are not allowed with an undirected network; use 'socialitycov'", call.=FALSE)

  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname
  if(is.numeric(x) && ! force.factor){
    # Numeric covariate.
    xm<-matrix(x,n,n,byrow=FALSE)
    cn<-paste("sendercov",attrname,sep=".")
    model <- .ergmm.add.fixed(model,xm,mean,var,cn)
  }else{
    # Factor covariate.
    x<-as.factor(x)
    ls<-levels(x)
    mean<-rep(mean,length.out=length(ls))
    var<-rep(var,length.out=length(ls))
    for(li in 2:length(ls)){
      l<-ls[[li]]
      xm<-matrix(x==l,n,n,byrow=FALSE)
      cn<-paste('sendercov',attrname,l,sep=".")
      model <- .ergmm.add.fixed(model, xm, mean[li], var[li], cn)
    }
  }
  model
}

#' @templateVar name receivercov
#' @title Receiver covariate effect
#' @description \emph{Deprecated for networks without self-loops. Use
#'   \code{\link{nodeicov-ergmTerm}},
#'   \code{\link{nodeifactor-ergmTerm}},
#'   \code{\link{nodecov-ergmTerm}} or
#'   \code{\link{nodefactor-ergmTerm}} instead.}
#'
#'   If the attribute is numeric, this term adds one covariate to the
#'   model equaling \code{attrname(i)}. If the attribute is not
#'   numeric or \code{force.factor==TRUE}, this term adds \eqn{p-1}
#'   covariates to the model, where \eqn{p} is the number of unique
#'   values of \code{attrname}.  The \eqn{k}th such covariate has the
#'   value \code{attrname(i) == value(k+1)}, where \code{value(k)} is
#'   the \eqn{k}th smallest unique value of the \code{attrname}
#'   attribute. This term only makes sense if the network is directed.
#'
#' @usage
#' # binary: receivercov(attrname, force.factor=FALSE, mean=0, var=9)
#'
#' # valued: receivercov(attrname, force.factor=FALSE, mean=0, var=9)
#'
#' @param attrname a character string giving the name of an attribute
#'   in the network's vertex attribute list.
#' @param force.factor logical, indicating if `attrname`'s value
#'   should be interpreted as categorical even if numeric.
#'
#' @template ergmTerm-latentnet-prior
#'
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#' @template ergmTerm-directed
#'
#' @concept dyad-independent
#' @concept directed
#' @concept latent
#' @concept loops
InitErgmTerm.receivercov <- nonlatent_error
InitErgmm.receivercov<-function (model, attrname, force.factor=FALSE, mean=0, var=9) 
{
  if(!has.loops(model[["Yg"]])) warning("Term `receivercov` is deprecated for networks without self-loops. Use `nodeicov`, `nodecov`, `nodeifactor`, or `nodefactor` from package `ergm` instead.")

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("receivercov() model term expected between 1 and 3 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  if (!is.directed(model[["Yg"]]))
    stop("Receiver covariates are not allowed with an undirected network; use 'socialitycov'", call.=FALSE)

  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname
  if(is.numeric(x) && ! force.factor){
    # Numeric covariate.
    xm<-matrix(x,n,n,byrow=TRUE)
    cn<-paste("receivercov",attrname,sep=".")
    model <- .ergmm.add.fixed(model,xm,mean,var,cn)
  }else{
    # Factor covariate.
    x<-as.factor(x)
    ls<-levels(x)
    mean<-rep(mean,length.out=length(ls))
    var<-rep(var,length.out=length(ls))
    for(li in 2:length(ls)){
      l<-ls[[li]]
      xm<-matrix(x==l,n,n,byrow=TRUE)
      cn<-paste('receivercov',attrname,l,sep=".")
      model <- .ergmm.add.fixed(model, xm, mean[li], var[li], cn)
    }
  }
  model
}

#' @templateVar name socialitycov
#' @title Sociality covariate effect
#' @description \emph{Deprecated for networks without self-loops. Use
#'   \code{\link{nodecov-ergmTerm}} or
#'   \code{\link{nodefactor-ergmTerm}} instead.}
#'
#'   If the attribute is numeric, this term adds one covariate to the
#'   model equaling \code{attrname(i)}. If the attribute is not
#'   numeric or \code{force.factor==TRUE}, this term adds \eqn{p-1}
#'   covariates to the model, where \eqn{p} is the number of unique
#'   values of \code{attrname}.  The \eqn{k}th such covariate has the
#'   value \code{attrname(i) == value(k+1)}, where \code{value(k)} is
#'   the \eqn{k}th smallest unique value of the \code{attrname}
#'   attribute. This term only makes sense if the network is directed.
#'
#' @usage
#' # binary: socialitycov(attrname, force.factor=FALSE, mean=0, var=9)
#'
#' # valued: socialitycov(attrname, force.factor=FALSE, mean=0, var=9)
#'
#' @param attrname a character string giving the name of an attribute
#'   in the network's vertex attribute list.
#' @param force.factor logical, indicating if `attrname`'s value
#'   should be interpreted as categorical even if numeric.
#'
#' @template ergmTerm-latentnet-prior
#'
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#'
#' @concept dyad-independent
#' @concept undirected
#' @concept bipartite
#' @concept latent
#' @concept loops
InitErgmTerm.socialitycov <- nonlatent_error
InitErgmm.socialitycov<-function (model, attrname, force.factor=FALSE, mean=0, var=9) 
{
  if(!has.loops(model[["Yg"]])) warning("Term `socialitycov` is deprecated for networks without self-loops. Use `nodeicov`, `nodecov`, `nodeifactor`, or `nodefactor` from package `ergm` instead.")

  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("socialitycov() model term expected between 1 and 4 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname
  if(is.numeric(x) && ! force.factor){
    # Numeric covariate.
    xm<-matrix(x,n,n,byrow=FALSE)+matrix(x,n,n,byrow=TRUE)
    cn<-paste("socialitycov",attrname,sep=".")
    model <- .ergmm.add.fixed(model,xm,mean,var,cn)
  }else{
    # Factor covariate.
    x<-as.factor(x)
    ls<-levels(x)
    mean<-rep(mean,length.out=length(ls))
    var<-rep(var,length.out=length(ls))
    for(li in 2:length(ls)){
      l<-ls[[li]]
      xm<-matrix(x==l,n,n,byrow=FALSE)+matrix(x==l,n,n,byrow=TRUE)
      cn<-paste('socialitycov',attrname,l,sep=".")
      model <- .ergmm.add.fixed(model, xm, mean[li], var[li], cn)
    }
  }
  model
}
