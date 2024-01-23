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
