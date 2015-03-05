\name{latentnet-package}
\alias{latentnet-package}
\alias{latentnet}
\title{Latent position and cluster models for networks}
\description{
  The package \code{latentnet} is used to fit latent cluster random effect
  models, where the probability of a network \eqn{g},
  on a set of nodes is a product of dyad probabilities, each of which is
  a GLM with linear component
  \eqn{\eta_{i,j}=\sum_{k=1}^p \beta_k X_{i,j,k}+d(Z_i,Z_j)+\delta_i+\gamma_j},
  where \eqn{X} is an array of dyad covariates, \eqn{\beta} is a vector
  of covariate coefficients, \eqn{Z_i} is the latent space position of
  node \eqn{i}, \eqn{d(\cdot,\cdot)} is a function of the two positions:
  either negative Euclidean (\eqn{-||Z_i-Z_j||}) or bilinear
  (\eqn{Z_i\cdot Z_j}), and \eqn{\delta} and \eqn{\gamma} are vectors of sender
  and receiver effects. (Note that these are different from the
  eigenmodel of Hoff (2007) ``Modeling homophily and stochastic
  equivalence in symmetric relational data'', fit by package \code{eigenmodel}.)

  The \code{\link{ergmm}} specifies models via: \code{g ~ <model terms>}
  where \code{g} is a \code{network} object For the list of possible
  \code{<model terms>}, see \link{terms.ergmm}. For the list of
  the possible dyad distribution families, see
  \link{families.ergmm}.
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
  Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum (2007).
 \emph{Model-Based Clustering for Social Networks}.
 Journal of the Royal Statistical Society: Series A (Statistics in
 Society), 170(2), 301-354.

 Peter D. Hoff (2005). \emph{Bilinear Mixed Effects Models for Dyadic
 Data}. Journal of the American Statistical Association, 100(469),
 286-295.

 Peter D. Hoff,  Adrian E. Raftery and Mark S. Handcock (2002).
 \emph{Latent space approaches to social network analysis}.
 Journal of the American Statistical Association, 97(460), 1090-1098.

 Pavel N. Krivitsky, Mark S. Handcock, Adrian E. Raftery, and Peter
 D. Hoff (2009).
 \emph{Representing degree distributions, clustering, and homophily in
   social networks with latent cluster random effects models}.
 Social Networks, 31(3), 204-213.

 Pavel N. Krivitsky and Mark S. Handcock (2008).
 \emph{Fitting Position Latent Cluster Models for Social Networks with
 \code{latentnet}}. Journal of Statistical Software, 24(5).

 Susan M. Shortreed, Mark S. Handcock, and Peter D. Hoff (2006).
 \emph{Positional Estimation within the Latent Space Model for
   Networks}. Methodology, 2(1), 24-33.
}
\keyword{graphs}
\keyword{package}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
\keyword{nonparametric}