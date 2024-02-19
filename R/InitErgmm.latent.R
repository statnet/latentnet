#  File R/InitErgmm.latent.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

InitErgmTerm.latent <- nonlatent_error
InitErgmm.latent <- function(model, d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
                             mean.var.mul=2, mean.var=NULL, pK.mul=1, pK=NULL){
  stop("Term ",sQuote("latent()")," has been deprecated in favor of ",sQuote("euclidean()"),".")
}

#' @templateVar name euclidean
#' @title Euclidean distance latent space, with optional clustering
#' @description Adds a term to the model equal to the negative
#'   Eucledean distance \eqn{-||Z_i-Z_j||}{-dist(Z[i],Z[j])}, where
#'   \eqn{Z_i}{Z[i]} and \eqn{Z_j}{Z[j]} are the positions of their
#'   respective actors in an unobserved social space. These positions
#'   may optionally have a finite spherical Gaussian mixture
#'   clustering structure. This term was previously called
#'   \code{latent}.
#'
#' @usage
#' # binary: euclidean(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
#' #             mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)
#'
#' # valued: euclidean(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
#' #             mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)
#'
#' @template ergmTerm-latentnet-latent-params
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#'
#' @concept dyad-independent
#' @concept undirected
#' @concept directed
#' @concept latent
InitErgmTerm.euclidean <- nonlatent_error
InitErgmm.euclidean<-function(model, d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
                           mean.var.mul=2, mean.var=NULL, pK.mul=1, pK=NULL){
  if(nargs()<2)
    stop(paste("euclidean() model term expected at least 1 argument, got ",
               nargs()-1, sep=""), call.=FALSE)
  if(d<=0) stop("Invalid latent space dimensionality given", call.=FALSE)
  if(!is.null(model[["d"]]) && model[["d"]]>0){
    stop("Only one latent position term can be added!")
  }

  model[["latent"]] <- "negative.Euclidean"
  
  model[["d"]] <- d
  model[["G"]] <- G

  model[["prior"]][["Z.var.mul"]]<-var.mul
  model[["prior"]][["Z.var"]]<-var
  model[["prior"]][["Z.var.df.mul"]]<-var.df.mul
  model[["prior"]][["Z.var.df"]]<-var.df
  model[["prior"]][["Z.mean.var.mul"]]<-mean.var.mul
  model[["prior"]][["Z.mean.var"]]<-mean.var
  model[["prior"]][["Z.pK.mul"]]<-pK.mul
  model[["prior"]][["Z.pK"]]<-pK
 
  if(!("Z.var" %in% names(model[["prior"]]))) model[["prior"]][["Z.var"]]<-model[["prior"]][["Z.var.mul"]]*(network.size(model[["Yg"]])/max(1,model[["G"]]))^(2/model[["d"]])
  if(!("Z.mean.var" %in% names(model[["prior"]]))) model[["prior"]][["Z.mean.var"]]<-model[["prior"]][["Z.mean.var.mul"]]*model[["prior"]][["Z.var"]]*max(1,model[["G"]])^(2/model[["d"]])
  if(!("Z.var.df" %in% names(model[["prior"]]))) model[["prior"]][["Z.var.df"]]<-model[["prior"]][["Z.var.df.mul"]]*sqrt(network.size(model[["Yg"]])/max(1,model[["G"]]))
  if(!("Z.pK" %in% names(model[["prior"]]))) model[["prior"]][["Z.pK"]]<-model[["prior"]][["Z.pK.mul"]]*sqrt(network.size(model[["Yg"]])/max(1,model[["G"]]))
  
  
  model
}

#' @templateVar name bilinear
#' @title Bilinear (inner-product) latent space, with optional clustering
#' @description Adds a term to the model equal to the inner product of
#'   the latent positions: \eqn{Z_i \cdot Z_j}{sum(Z[i]*Z[j])}, where
#'   \eqn{Z_i}{Z[i]} and \eqn{Z_j}{Z[j]} are the positions of their
#'   respective actors in an unobserved social space. These positions
#'   may optionally have a finite spherical Gaussian mixture
#'   clustering structure. \emph{Note: For a bilinear latent space
#'   effect, two actors being closer in the clustering sense does not
#'   necessarily mean that the expected value of a tie between them is
#'   higher. Thus, a warning is printed when this model is combined
#'   with clustering.}
#'
#' @usage
#' # binary: bilinear(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
#' #             mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)
#'
#' # valued: bilinear(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
#' #             mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)
#'
#' @template ergmTerm-latentnet-latent-params
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#'
#' @concept dyad-independent
#' @concept undirected
#' @concept directed
#' @concept latent
InitErgmTerm.bilinear <- nonlatent_error
InitErgmm.bilinear<-function(model, d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
                           mean.var.mul=2, mean.var=NULL, pK.mul=1, pK=NULL){
  if(nargs()<2)
    stop(paste("bilinear() model term expected at least 1 argument, got ",
               nargs()-1, sep=""), call.=FALSE)
  if(d<=0) stop("Invalid latent space dimensionality given", call.=FALSE)
  if(!is.null(model[["d"]]) && model[["d"]]>0){
    stop("Only one latent position term can be added!")
  }

  model[["latent"]] <- "bilinear"
  
  model[["d"]] <- d
  model[["G"]] <- G

  if(G>0) warning("For a bilinear (inner-product) latent position effect, two actors being closer in a clustering sense does not necessarily mean a higher expected value of their relationship. Thus, clustering might not be interpretable.")
  
  model[["prior"]][["Z.var.mul"]]<-var.mul
  model[["prior"]][["Z.var"]]<-var
  model[["prior"]][["Z.var.df.mul"]]<-var.df.mul
  model[["prior"]][["Z.var.df"]]<-var.df
  model[["prior"]][["Z.mean.var.mul"]]<-mean.var.mul
  model[["prior"]][["Z.mean.var"]]<-mean.var
  model[["prior"]][["Z.pK.mul"]]<-pK.mul
  model[["prior"]][["Z.pK"]]<-pK
  
  if(!("Z.var" %in% names(model[["prior"]]))) model[["prior"]][["Z.var"]]<-model[["prior"]][["Z.var.mul"]]*(network.size(model[["Yg"]])/max(1,model[["G"]]))^(2/model[["d"]])
  if(!("Z.mean.var" %in% names(model[["prior"]]))) model[["prior"]][["Z.mean.var"]]<-model[["prior"]][["Z.mean.var.mul"]]*model[["prior"]][["Z.var"]]*max(1,model[["G"]])^(2/model[["d"]])
  if(!("Z.var.df" %in% names(model[["prior"]]))) model[["prior"]][["Z.var.df"]]<-model[["prior"]][["Z.var.df.mul"]]*sqrt(network.size(model[["Yg"]])/max(1,model[["G"]]))
  if(!("Z.pK" %in% names(model[["prior"]]))) model[["prior"]][["Z.pK"]]<-model[["prior"]][["Z.pK.mul"]]*sqrt(network.size(model[["Yg"]])/max(1,model[["G"]]))  
  
  model
}

#' @templateVar name euclidean2
#' @title Squared euclidean distance latent space, with optional clustering
#' @description Adds a term to the model equal to the negative
#'   Eucledean distance \eqn{-||Z_i-Z_j||^2}{-dist(Z[i],Z[j])^2}, where
#'   \eqn{Z_i}{Z[i]} and \eqn{Z_j}{Z[j]} are the positions of their
#'   respective actors in an unobserved social space. These positions
#'   may optionally have a finite spherical Gaussian mixture
#'   clustering structure. This term was previously called
#'   \code{latent}.
#'
#' @usage
#' # binary: euclidean(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
#' #             mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)
#'
#' # valued: euclidean(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
#' #             mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)
#'
#' @template ergmTerm-latentnet-latent-params
#' @template ergmTerm-general
#' @template ergmTerm-latentnet-general
#'
#' @concept dyad-independent
#' @concept undirected
#' @concept directed
#' @concept latent
InitErgmTerm.euclidean2 <- nonlatent_error
InitErgmm.euclidean2<-function(model, d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
                           mean.var.mul=2, mean.var=NULL, pK.mul=1, pK=NULL){
  if(nargs()<2)
    stop(paste("euclidean2() model term expected at least 1 argument, got ",
               nargs()-1, sep=""), call.=FALSE)
  if(d<=0) stop("Invalid latent space dimensionality given", call.=FALSE)
  if(!is.null(model[["d"]]) && model[["d"]]>0){
    stop("Only one latent position term can be added!")
  }

  model[["latent"]] <- "negative.Euclidean2"
  
  model[["d"]] <- d
  model[["G"]] <- G

  model[["prior"]][["Z.var.mul"]]<-var.mul
  model[["prior"]][["Z.var"]]<-var
  model[["prior"]][["Z.var.df.mul"]]<-var.df.mul
  model[["prior"]][["Z.var.df"]]<-var.df
  model[["prior"]][["Z.mean.var.mul"]]<-mean.var.mul
  model[["prior"]][["Z.mean.var"]]<-mean.var
  model[["prior"]][["Z.pK.mul"]]<-pK.mul
  model[["prior"]][["Z.pK"]]<-pK
 
  if(!("Z.var" %in% names(model[["prior"]]))) model[["prior"]][["Z.var"]]<-model[["prior"]][["Z.var.mul"]]*(network.size(model[["Yg"]])/max(1,model[["G"]]))^(2/model[["d"]])
  if(!("Z.mean.var" %in% names(model[["prior"]]))) model[["prior"]][["Z.mean.var"]]<-model[["prior"]][["Z.mean.var.mul"]]*model[["prior"]][["Z.var"]]*max(1,model[["G"]])^(2/model[["d"]])
  if(!("Z.var.df" %in% names(model[["prior"]]))) model[["prior"]][["Z.var.df"]]<-model[["prior"]][["Z.var.df.mul"]]*sqrt(network.size(model[["Yg"]])/max(1,model[["G"]]))
  if(!("Z.pK" %in% names(model[["prior"]]))) model[["prior"]][["Z.pK"]]<-model[["prior"]][["Z.pK.mul"]]*sqrt(network.size(model[["Yg"]])/max(1,model[["G"]]))
  
  
  model
}
