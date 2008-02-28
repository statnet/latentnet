\name{latentnet}
\alias{latentnet}
\title{Latent position and cluster models for networks}
\description{
  The package \code{latentnet} is used to fit latent cluster 
  models, where the probability of a network \eqn{g},
  on a set of nodes is a product of dyad probabilities, each of which is
  a GLM with linear component
  \eqn{\eta_{i,j}=\sum_{k=1}^p \beta_k X_{i,j,k}+||Z_i-Z_j||},
  where \eqn{X} is an array of dyad covariates, \eqn{\beta} is a vector
  of covariate coefficients and \eqn{Z_i} is the latent space position of
  node \eqn{i}.

  The \code{\link{ergmm}} specifies models via: \code{g ~ <model terms>}
  where \code{g} is a \code{network} object
  For the list of possible \code{<model terms>}, see 
  \code{\link{terms.ergmm}}.
}

\details{

    The arguments in the \code{\link{ergmm}} function specific to latent 
    variable models
    are \code{ergmm.control}. See the help page for \code{\link{ergmm}} for the 
    details.

    The result of a latent variable model fit is an \code{\link{ergmm}} object.
    Hence the \code{\link{summary}}, \code{print}, and \code{plot} functions
    apply to the fits.
    The \code{\link{plot.ergmm}} function has many options specific to latent 
    variable models.
% are \code{mle}, \code{pie}, \code{contour},
% \code{contour.colors}.
    See the help page for \code{\link{plot.ergmm}} for the details.
}
\value{\code{\link{ergmm}} returns an object of class 'ergmm' that is a list.
}
\seealso{\code{\link{ergmm}}, \link{terms.ergmm}}
\references{
 Peter D. Hoff,  Adrian E. Raftery and Mark S. Handcock.
 \emph{Latent space approaches to social network analysis}.
 Journal of the American Statistical Association, Dec 2002, Vol. 97, 
 Iss. 460;  pg. 1090-1098.

 Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum.
 \emph{Model-Based Clustering for Social Networks}.
 Journal of the Royal Statistical Society: Series A (Statistics in
 Society), 2007, Vol. 170, Iss. 2; pg. 301-354.

 Susan M. Shortreed, Mark S. Handcock, and Peter D. Hoff.
 \emph{Positional Estimation within the Latent Space Model for
   Networks}. Methodology, 2006, Vol. 2, Iss. 1; pg. 24-33.
}
\keyword{graphs}
\keyword{package}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
\keyword{nonparametric}