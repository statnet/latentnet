%  File man/ergmm.par.list.object.Rd in package latentnet, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergmm.par.list.object}
\alias{ergmm.par.list.object}
\alias{ergmm.par.list}
\alias{as.ergmm.par.list}
\alias{as.mcmc.list.ergmm.par.list}
\alias{[.ergmm.par.list}
\alias{[[.ergmm.par.list}
\alias{$.ergmm.par.list}
\alias{length.ergmm.par.list}
\alias{stack.ergmm.par.list.list}
\alias{unstack.ergmm.par.list}
\alias{del.iteration}
\title{A List of ERGMM Parameter Configuration}
\description{
  A class \code{\link[=ergmm.par.list.object]{ergmm.par.list}} to represent a
  series of parameter configurations for the same exponential random graph
  mixed model.
}
\details{
  \code{[[} operator with a
  numeric or integer index returns a
  list with the the
  configuration with that index. \code{[} operator given a numeric
  vector returns a \code{ergmm.par.list} object with the subset of
  configurations with the indices given.
  
  The structure of \code{ergmm.par.list} is derived from named lists, with each entry having an
  additional dimension (always the first one), indexed by
  configuration. That is, scalars become vectors, vectors become
  matrixes with the original vectors in rows, and matrixes become
  3-dimensional arrays, with the original matrixes indexed by their
  first dimension. See \code{\link{terms.ergmm}} for comon elements of
  these configurations.

  In some cases, such as when representing MCMC or optimization output,
  the object may also have some of the following elements:
  \describe{
    \item{\code{mlp}}{\eqn{\log
	p(Y,Z,\beta,\mu,\sigma,\delta,\gamma,\sigma_\delta,\sigma_\gamma,|K)},
      the joint
      probability/density of network, the covariate coefficients, the
      latent space positions and parameters, and the random effects and
      their variances, conditional on cluster assignments.}
    \item{\code{lpY}}{\eqn{\log p(Y|\dots)}, depending on the model, the log-probability or
      log-density of the network conditional on all the parameters.}
    \item{\code{lpZ}}{\eqn{\log p(Z|\mu,\sigma,K)}, the log-density of latent space positions conditional on
      latent space or cluster parameters and cluster assignments.}
    \item{\code{lpbeta}}{\eqn{\log p(\beta)}, the prior log-density of
      the covariate coefficients.}
    \item{\code{lpRE}}{\eqn{\log
	p(\delta,\gamma|\sigma_\delta,\sigma_\gamma)}, the log-density of all random effects, conditional on
      their respective variances.}
    \item{\code{lpLV}}{\eqn{\log p(\mu,\sigma)}, the prior log-density
      of latent space or cluster parameters (but not that of the cluster
      assignments).}
    \item{\code{lpREV}}{\eqn{\log p(\sigma_\delta,\sigma_\gamma)}, the prior log-density of all random effect variances.}
    \item{\code{Z.rate}}{Proportion of single-vertex proposals accepted over the preceding
      interval.}
    \item{\code{beta.rate}}{
      Proportion of group proposals accepted over the preceding interval.}
  }
}
\seealso{
  \code{\link{ergmm}}
}
\keyword{graphs}
\keyword{utilities}
\keyword{methods}
\keyword{list}
\keyword{manip}