\name{ergmm.control}
\alias{ergmm.control}
\title{Auxiliary for Controlling ERGMM Fitting}
\description{
  Auxiliary function as user interface for \code{ergmm}
  fitting. Typically only used when calling \code{ergmm}. It is used to
  set parameters that affect the sampling but do not affect the posterior distribution.
}
\usage{
ergmm.control(sample.size=4000,
              burnin=10000,
              interval=10,
              threads=1,
              mle.maxit=100,
              Z.delta=0.6,
              RE.delta=0.6,
              group.deltas=0.4,
              pilot.runs=4,
              pilot.factor=0.8,
              pilot.discard.first=0.5,
              target.acc.rate=0.234,
              backoff.threshold=0.05,
              backoff.factor=0.1,
              accept.all=FALSE,
              store.burnin=FALSE)
}
\arguments{
  \item{sample.size}{The number of draws to be taken from the posterior distribution.}
  \item{burnin}{The number of initial MCMC iterations to be discarded.}
  \item{interval}{The number of iterations between consecutive draws.}
  \item{threads}{The number of chains to run. If greater than 1,
    package \code{\link[snowFT:snowFT-package]{snowFT}} is used to take
    advantage of any multiprocessing or distributed computing
    capabilities that may be available. Currently, only PVM (via
    \code{rpvm}) has been tested. Note, also, that PVM daemon needs to
    be started before the package is loaded. } 
  \item{mle.maxit}{Maximum number of iterations for computing the
    starting values, posterior modes, MLEs, MKL estimates, etc..}
  \item{Z.delta}{Standard deviation of the proposal for the jump in the
    individual latent space position, or its starting value for the tuner.}
  \item{RE.delta}{Standard deviation of the proposal for the jump in the
    individual random effects values, or its starting value for the
    tuner.}
  \item{group.deltas}{A scalar, a vector, or a matrix of an
    appropriate size, giving the initial proposal structure for the
    ``group proposal'' of a jump in covariate coefficients, scaling of
    latent space positions, and a shift in random ffects. If a matrix of
    an appropriate size is given, it is used as a matrix of coefficients
    for a correlated proposal. If a vector is given, an independent
    proposal is used with the corresponding elements being proposal
    standard deviations. If a scalar is given, it is used as a
    multiplier for an initial heuristic for the proposal structure. It
    is usually best to leave this argument alone and let the adaptive
    sampling be used.}
  \item{pilot.runs}{Number of pilot runs into which to split the
    burn-in period. After each pilot run, the proposal standard
    deviations and coefficients \code{Z.delta}, \code{RE.delta}, and
    \code{group.deltas} are reevaluated. If set to \code{0}, disables
    adaptive sampling, and only makes a single burn-in run.}
  \item{pilot.factor}{Initial value for the factor by which the
    coefficients gotten by a Choletsky decomposition of the pilot sample
    covariance matrix are multiplied.}
  \item{pilot.discard.first}{Proportion of draws from the pilot run to
    discard for estimating acceptance rate and group proposal
    covariance.}
  \item{target.acc.rate}{Taget acceptance rate for the proposals
    used. After a pilot run, the proposal variances are adjusted upward
    if the acceptance rate is above this, and downward if below.}
  \item{backoff.threshold}{If a pilot run's acceptance rate is below
    this, redo it with drastically reduced proposal standard deviation.
    Set to \code{0} to disable this behavior.}
  \item{backoff.factor}{Factor by which to multiply the relevant
    proposal standard deviation if its acceptance rate fell below the
    backoff threshold. }
  \item{accept.all}{Forces all proposals to be accepted
    unconditionally. Use only in debugging proposal distributions!}
  \item{store.burnin}{If \code{TRUE}, the samples from the burnin are
    also stored and returned, to be used in MCMC diagnostics.}
}
\value{
  A list with the arguments as components.
}
\seealso{ \code{\link{ergmm}} }
\examples{
data(sampson)
## Shorter run than default.
samp.fit<-ergmm(samplike~latent(d=2,G=3)+rreceiver,
control=ergmm.control(burnin=1000,sample.size= 2000,interval=5))
}
\keyword{graphs}
\keyword{misc}