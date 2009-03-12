\name{ergmm.object}
\alias{ergmm.object}
\title{Class of Fitted Exponential Random Graph Mixed Models}
\description{
  A class \code{\link[=ergmm.object]{ergmm}} to represent a fitted exponential random graph
  mixed model. The output of \code{\link{ergmm}}.
}
\details{
  There are methods \code{\link{summary.ergmm}}, \code{print.ergmm},
  \code{\link{plot.ergmm}}, \code{\link{predict.ergmm}},
  \code{\link{as.mcmc.ergmm}}, and \code{\link{as.mcmc.list.ergmm}}.
  
  The structure of \code{ergmm} is as follows:
  \tabular{ll}{
    \code{sample}\tab An object of class
    \code{\link[=ergmm.par.list.object]{ergmm.par.list}} containing the
    MCMC sample from the posterior. If the run had multiple threads, their output is concatenated.\cr
    \code{mcmc.mle}\tab A list containing the parameter
    configuration of the highest-likelihood MCMC iteration. \cr
    \code{mcmc.pmode}\tab A list containing the parameter
    configuration of the highest-joint-density (conditional on cluster
    assignments) MCMC iteration. \cr
    \code{mkl}\tab A list containing the MKL estimate. \cr
    \code{model}\tab A list containing the model
    that was fitted.\cr
    \code{prior}\tab A list containing the
    information about the prior distribution used. It can be passed as
    parameter \code{prior} to \code{\link{ergmm}} to reproduce the prior
    in a new fit.\cr
    \code{control}\tab A list containing the
    information about the model fit settings that do not affect the
    posterior distribution. It can be passed as
    parameter \code{control} to \code{\link{ergmm}} to reproduce control
    parameters in a new fit.\cr
    \code{mle}\tab A list containing the MLE,
    conditioned on cluster assignments.\cr
    \code{pmode}\tab A list containing the posterior mode,
    conditioned on cluster assignments.\cr
    \code{burnin.start}\tab A list containing the starting
    value for the burnin.\cr
    \code{main.start}\tab  A list (or a list of lists, for a
    multithreaded run) containing the starting
    value for the sampling. \cr
  }
}
\seealso{
  \code{\link{ergmm}},
  \code{\link{summary.ergmm}},
  \code{\link{plot.ergmm}},
  \code{\link{predict.ergmm}},
  \code{\link{as.mcmc.ergmm}},
  \code{\link{as.mcmc.list.ergmm}}
}
\keyword{graphs}
\keyword{regression}
\keyword{models}