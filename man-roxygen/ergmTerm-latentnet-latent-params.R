#'
#' @param d The dimension of the latent space.
#' @param G The number of groups (0 for no clustering).
#' @param var.mul In the absence of \code{var}, this argument will be
#'   used as a scaling factor for a function of average cluster size
#'   and latent space dimension to set \code{var}. To set it in the
#'   \code{prior} argument to \code{\link{ergmm}}, use
#'   \code{Z.var.mul}.
#' @param var If given, the scale parameter for the
#'   scale-inverse-chi-squared prior distribution of the
#'   within-cluster variance. To set it in the \code{prior} argument
#'   to \code{\link{ergmm}}, use \code{Z.var}.
#' @param var.df.mul In the absence of \code{var.df}, this argument is
#'   the multiplier for the square root of average cluster size, which
#'   serves in place of \code{var.df}. To set it in the \code{prior}
#'   argument to \code{\link{ergmm}}, use \code{Z.var.df.mul}.
#' @param var.df The degrees of freedom parameter for the
#'   scale-inverse-chi-squared prior distribution of the
#'   within-cluster variance. To set it in the \code{prior} argument
#'   to \code{\link{ergmm}}, use \code{Z.var.df}.
#' @param mean.var.mul In the absence of \code{mean.var}, the
#'   multiplier for a function of number of vertices and latent space
#'   dimension to set \code{mean.var}. To set it in the \code{prior}
#'   argument to \code{\link{ergmm}}, use \code{Z.mean.var.mul}.
#' @param mean.var The variance of the spherical Gaussian prior
#'   distribution of the cluster means. To set it in the \code{prior}
#'   argument to \code{\link{ergmm}}, use \code{Z.mean.var}.
#' @param pK.mul In the absence of \code{pK}, this argument is the
#'   multiplier for the square root of the average cluster size, which
#'   is used as \code{pK}. To set it in the \code{prior} argument to
#'   \code{\link{ergmm}}, use \code{Z.pK}.
#' @param pK The parameter of the Dirichilet prior distribution of
#'   cluster assignment probabilities. To set it in the \code{prior}
#'   argument to \code{\link{ergmm}}, use \code{Z.pK}.
#'
#' @details The following parameters are associated with this term:
#'   \describe{
#'   \item{\code{Z}}{ Numeric matrix with rows being latent space
#'   positions.}
#'   \item{\code{Z.K} (when \eqn{\code{G}>0})}{ Integer vector of
#'   cluster assignments. }
#'   \item{\code{Z.mean} (when \eqn{\code{G}>0})}{ Numeric matrix
#'   with rows being cluster means. }
#'   \item{\code{Z.var} (when \eqn{\code{G}>0})}{ Depending on the
#'   model, either a numeric vector with within-cluster variances
#'   or a numeric scalar with the overal latent space variance. }
#'   \item{\code{Z.pK} (when \eqn{\code{G}>0})}{ Numeric vector of
#'   probabilities of a vertex being in a particular cluster.}
#'   }
