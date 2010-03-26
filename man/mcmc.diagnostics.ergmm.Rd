\name{mcmc.diagnostics.ergmm}
\alias{mcmc.diagnostics.ergmm}
\alias{mcmc.diagnostics}
\alias{mcmc.diagnostics.default} 
\title{Conduct MCMC diagnostics on an ERGMM fit}
\description{
  This function creates simple diagnostic plots for the MCMC sampled
  statistics produced from a fit. It also prints the Raftery-Lewis 
  diagnostics, indicates if they are sufficient,
  and suggests the run length required.
}
\usage{
\method{mcmc.diagnostics}{ergmm}(x,which.diags=c("cor","acf","trace","raftery"),
                                 burnin=FALSE,
                                 which.vars=NULL,
                                 vertex.i=c(1),\dots)
}
\arguments{
  \item{x}{ An object of class \code{\link[=ergmm.object]{ergmm}}. }
  \item{which.diags}{A list of diagnostics to produce. "cor" is the
    correlation matrix of the statistics, "acf" plots the
    autocorrelation functions, "trace" produces trace plots and density
    estimates, and "raftery" produces the Raftery-Lewis statistics.}
  \item{burnin}{ If not \code{FALSE}, rather than perform diagnostics on
    the sampling run, performs them on the pilot run whose index is given. }
  \item{which.vars}{A named list mapping variable names to the indices
    to include. If given, overrides the defaults and all arguments that follow.}
  \item{vertex.i}{ A numeric vector of vertices whose latent space
    coordinates and random effects to include. }
  \item{\dots}{Additional arguments. None are supported at the moment.}
}

\details{
  Produces the plots per \code{which.diags}.
  Autocorrelation function that is printed if "acf" is requested is for
  lags \code{0} and \code{interval}.
}
\value{
  \code{mcmc.diagnostics.ergmm} returns a table of Raftery-Lewis diagnostics.
}

\seealso{\code{\link{ergmm}}, \code{\link{ergmm.object}},
  \code{\link[coda]{raftery.diag}},
  \code{\link[coda]{autocorr}},
  \code{\link[coda]{plot.mcmc.list}}
}
\examples{
#
data(sampson)
#
# test the mcmc.diagnostics function
#
gest <- ergmm(samplike ~ euclidean(d=2),
              control=ergmm.control(burnin=1000,interval=5))
summary(gest)
#
# Plot the traces and densities
#
mcmc.diagnostics(gest)
}
\keyword{graphs}
\keyword{hplot}
\keyword{debugging}