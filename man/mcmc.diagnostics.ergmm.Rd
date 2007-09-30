\name{mcmc.diagnostics.ergmm}
\alias{mcmc.diagnostics.ergmm}
\alias{mcmc.diagnostics}
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
                                 vertex.i=c(1))
}
\arguments{
  \item{x}{ An object of class \code{\link[=ergmm.object]{ergmm}}. }
  \item{which.diags}{A list of diagnostics to produce. "cor" is the correlation
matrix of the statistics, "acf" plots the autocorrelation functions, "trace"
produces trace plots, and "raftery" produces the Raftery-Lewis statistics.} 
  \item{burnin}{ If \code{TRUE}, generates an
    \code{\link[coda]{mcmc.list}} object for the burnin (if stored)
    instead of the main sampling run. }
  \item{which.vars}{A named list mapping variable names to the indices
    to include. If given, overrides the defaults and all arguments that follow.}
  \item{vertex.i}{ A numeric vector of vertices whose latent space
    coordinates and random effects to include. }
}

\details{
The plots produced are a trace of the sampled output
and a density estimate for each variable in the chain.

Autocorrelation with lags 0 and 1 and a Raftery-Lewis diagnostic is prited.
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
gest <- ergmm(samplike ~ latent(d=2))
summary(gest)
#
# Plot the traces and densities
#
mcmc.diagnostics(gest)
}
\keyword{graphs}
