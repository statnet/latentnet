#  File R/InitErgmm.random.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

#' @templateVar name rsender
#' @templateVar reffect sender
#' @title Random sender effect
#' @description Adds a random sender effect to the model, with normal
#'   prior centered around \eqn{0}{0} and a variance that is
#'   estimated. Can only be used on directed networks.
#'
#' @usage
#' # binary: rsender(var=1, var.df=3)
#'
#' # valued: rsender(var=1, var.df=3)
#'
#' @template ergmTerm-latentnet-random
#' @template ergmTerm-general
#' @template ergmTerm-directed
#'
#' @concept dyad-independent
#' @concept directed
#' @concept latent
InitErgmTerm.rsender <- nonlatent_error
InitErgmm.rsender<-function(model, var=1, var.df=3){
  if (!is.directed(model[["Yg"]]))
    stop("Sender effects are not allowed with an undirected network; use 'sociality'", call.=FALSE)
  model[["sender"]]<-TRUE
  model[["prior"]][["sender.var"]]<-var
  model[["prior"]][["sender.var.df"]]<-var.df
  model
}

#' @templateVar name rreceiver
#' @templateVar reffect receiver
#' @title Random receiver effect
#' @description Adds a random receiver effect to the model, with normal
#'   prior centered around \eqn{0}{0} and a variance that is
#'   estimated. Can only be used on directed networks.
#'
#' @usage
#' # binary: rreceiver(var=1, var.df=3)
#'
#' # valued: rreceiver(var=1, var.df=3)
#'
#' @template ergmTerm-latentnet-random
#' @template ergmTerm-general
#' @template ergmTerm-directed
#'
#' @concept dyad-independent
#' @concept directed
#' @concept latent
InitErgmTerm.rreceiver <- nonlatent_error
InitErgmm.rreceiver<-function(model, var=1, var.df=3){
  if (!is.directed(model[["Yg"]]))
    stop("receiver effects are not allowed with an undirected network; use 'sociality'", call.=FALSE)
  model[["receiver"]]<-TRUE
  model[["prior"]][["receiver.var"]]<-var
  model[["prior"]][["receiver.var.df"]]<-var.df
  model
}

#' @templateVar name rsociality
#' @templateVar reffect sociality
#' @title Random sociality effect
#' @description Adds a random sociality effect to the model, with normal
#'   prior centered around \eqn{0}{0} and a variance that is
#'   estimated. Can only be used on directed networks.
#'
#' @usage
#' # binary: rsociality(var=1, var.df=3)
#'
#' # valued: rsociality(var=1, var.df=3)
#'
#' @template ergmTerm-latentnet-random
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept undirected
#' @concept directed
#' @concept latent
InitErgmTerm.rsociality <- nonlatent_error
InitErgmm.rsociality<-function(model, var=1, var.df=3){
  model[["sociality"]]<-TRUE
  model[["prior"]][["sociality.var"]]<-var
  model[["prior"]][["sociality.var.df"]]<-var.df
  model
}
