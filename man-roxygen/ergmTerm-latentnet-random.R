#  File man-roxygen/ergmTerm-latentnet-random.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#'
#' @param var The scale parameter for the scale-inverse-chi-squared
#'   prior distribution of the <%= reffect %> effect variance. To set
#'   it in the \code{prior} argument to \code{\link{ergmm}}, use
#'   \code{<%= reffect %>.var}.
#'
#' @param var.df The degrees of freedom parameter for the
#'   scale-inverse-chi-squared prior distribution of the <%= reffect %> effect
#'   variance. To set it in the \code{prior} argument to
#'   \code{\link{ergmm}}, use \code{<%= reffect %>.var.df}.
#'
#'
#' @details The following parameters are associated with this term:
#'   \describe{
#'
#'   \item{\code{<%= reffect %>}}{ Numeric vector of values of each
#' 	  vertex's random <%= reffect %> effect.}
#'
#'   \item{\code{<%= reffect %>.var}}{ Random <%= reffect %> effect's variance.}
#'
#' }
