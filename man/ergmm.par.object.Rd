\name{ergmm.par.object}
\alias{ergmm.par}
\alias{ergmm.par.object}
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
  \tabular{ll}{
    \code{beta}\tab Numeric vector of covariate coefficients. \cr
    \code{Z.K}\tab Integer vector of cluster assignments. \cr
    \code{Z.mean}\tab Numeric matrix with rows being cluster means. \cr
    \code{Z.var}\tab Depending on the model, either a numeric vector with
    within-cluster variances or a numeric scalar with the overal latent space variance. \cr
    \code{Z.pK}\tab Numeric vector of probabilities of a vertex being
    in a particular cluster. \cr
    \code{Z}\tab Numeric matrix with rows being latent space positions. \cr
    \code{sender, receiver, sociality}\tab Numeric vector of values of each
    vertex's respective random effect.\cr
    \code{sender.var, receiver.var, sociality.var}\tab Numeric scalar
    with the respective random effect's variance.\cr
  }
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