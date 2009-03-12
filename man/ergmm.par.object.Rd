\name{ergmm.par.object}
\alias{ergmm.par}
\alias{ergmm.par.object}
\alias{as.ergmm.par}
\alias{[[.ergmm.par}
\alias{$.ergmm.par}
\alias{ergmm.par}
\title{An ERGMM Parameter Configuration}
\description{
  A class \code{\link[=ergmm.par.object]{ergmm.par}} to represent a
  parameter configuration for an exponential random graph
  mixed model.
}
\details{
  An \code{ergmm.par} object is essentially a named list of parameter
  values. It can be accessed in much the same way a list can, but no
  partial matching is performed.
  
  The structure of \code{ergmm.par} is a list which may include some of the following:

  In some cases (such as when representing MCMC or optimization output),
  the object may also have some of the following:
  \tabular{ll}{
    \code{mlp}: \eqn{\log
      p(Y,Z,\beta,\mu,\sigma,\delta,\gamma,\sigma_\delta,\sigma_\gamma,|\dots)}\tab Joint
    probability/density of network, the covariate coefficients, the
    latent space positions and parameters, and the random effects and
    their variances, conditional on cluster assignments. \cr
    \code{llk}: \eqn{\log p(Y|\dots)}\tab Depending on the model, the log-probability or
    log-density of the network conditional on all the parameters. \cr
    \code{lpZ}: \eqn{\log p(Z|\mu,\sigma,K)} \tab Log-density of latent space positions conditional on
    latent space or cluster parameters and cluster assignments. \cr
    \code{lpbeta}: \eqn{\log p(\beta)}\tab Prior log-density of the covariate coefficients. \cr
    \code{lpRE}: \eqn{\log p(\delta,\gamma|\sigma_\delta,\sigma_\gamma)}\tab Log-density of all random effects, conditional on
    their respective variances. \cr
    \code{lpLV}: \eqn{\log p(\mu,\sigma)}\tab Prior log-density of latent space or cluster
    parameters (but not that of the cluster assignments).\cr
    \code{lpREV}: \eqn{\log p(\sigma_\delta,\sigma_\gamma)}\tab
    Prior log-density of all random effect variances. \cr
    \code{Z.rate}\tab
    Proportion of single-vertex proposals accepted over the preceding
    interval. \cr
    \code{beta.rate}\tab
    Proportion of group proposals accepted over the preceding interval. \cr
  }
}
\seealso{
  \code{\link{ergmm.par.list}}
}
\keyword{graphs}
\keyword{utilities}
\keyword{methods}
\keyword{manip}