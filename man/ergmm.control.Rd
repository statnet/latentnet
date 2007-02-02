\name{ergmm.control}
\alias{ergmm.control}
\title{Auxiliary for Controlling ERGMM Fitting}
\description{
  Auxiliary function as user interface for \code{ergmm}
  fitting. Typically only used when calling \code{ergmm}. It is used to
  set parameters that affect the sampling but do not affect the posterior distribution.
}
\usage{
ergmm.control(samplesize = 5000, burnin = 1000, interval = 10, threads = 1, mle.maxit = 400, tune = FALSE, tuning.runs = 40, tuning.runsize = 4, Z.delta = 0.4, beta.delta = 0.4, Z.scl.delta = 0.02, RE.delta = 0.3, flyapart.penalty = 0)
}
\arguments{
  \item{samplesize}{The number of draws to be taken from the posterior distribution.}
  \item{burnin}{The number of MCMC iterations to be discarded.}
  \item{interval}{The number of iterations between consecutive draws.}
  \item{threads}{The number of chains to run. If greater than 1,
    \code{\link{snowFT}} is used to take advantage of any
    multiprocessing or distributed computing capabilities that may be available.}
  \item{mle.maxit}{Maximum number of iterations for computing the
    starting values, posterior modes, MLEs, MKL estimates, etc..}
  \item{tune}{If \code{TRUE}, the proposal variances are tuned before
    burnin and again before sampling begins. The tuner is somewhat experimental.}
  \item{tuning.runs}{Number of points at which to evaluate tuning.}
  \item{tuning.runsize}{Number of iterations per point.}
  \item{Z.delta}{Standard deviation of the proposal for the jump in the
    individual latent space position, or its starting value for the tuner.}
  \item{beta.delta}{Standard deviation of the proposal for the jump in the
    covariate coefficients, or its starting value for the tuner.}
  \item{Z.scl.delta}{Standard deviation of the log of the proposal for
    the scaling of latent space positions, or its starting value for the tuner.}
  \item{RE.delta}{Standard deviation of the proposal for the jump in the
    individual random effects values, or its starting value for the tuner.}
  \item{flyapart.penalty}{Multiplier for \eqn{\sum{\theta_i}^2}, that
    is, the sum of squares of all the parameter values, used to
    penaltize the log-posterior-density when finding the MCMC starting values.}
}
\value{
  A list with the arguments as components.
}
\author{Pavel N. Krivitsky \email{pavel@stat.washington.edu}}
\seealso{ \code{\link{ergmm}} }
\examples{
data(sampson)
## Longer run than default.
ergmm(samplike~latent(d=2,G=3)+rreceiver,control=ergmm.control(samplesize = 20000))
}
\keyword{ graphs }

