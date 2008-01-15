\name{as.mcmc.list.ergmm}
\alias{as.mcmc.list.ergmm}
\alias{as.mcmc.ergmm}
\title{Convert an ERGMM Object to an MCMC list object for Diagnostics.}
\description{
  A function to extract a subset of MCMC samples from an object of class
  \code{\link[=ergmm.object]{ergmm}} to an \code{\link[coda]{mcmc.list}} object.
}
\usage{
\method{as.mcmc}{ergmm}(x, burnin = FALSE, which.vars = NULL, vertex.i = c(1))
\method{as.mcmc.list}{ergmm}(x, burnin = FALSE, which.vars = NULL, vertex.i = c(1))
}
\arguments{
  \item{x}{ An object of class \code{\link[=ergmm.object]{ergmm}}. }
  \item{burnin}{ If \code{TRUE}, generates an
    \code{\link[coda]{mcmc.list}} object for the burnin (if stored)
    instead of the main sampling run. }
  \item{which.vars}{A named list mapping variable names to the indices
    to include. If given, overrides the defaults and all arguments that follow.}
  \item{vertex.i}{ A numeric vector of vertices whose latent space
    coordinates and random effects to include. }
}
\details{
  Unless \code{which.vars} is specified, the
  \code{\link[coda]{mcmc.list}} returned also includes all of the
  covariate coefficients.
  
  Regardless of whether the MCMC run was single- or multi-threaded, this
  function returns an \code{\link[coda]{mcmc.list}}, with a single
  thread, if necessary.
}
\value{
  A \code{\link[coda]{mcmc.list}} object with the samples of the
  selected subset of the variables.
}
\author{ Pavel N. Krivitsky \email{pavel@stat.washington.edu}}
\seealso{ \code{\link{ergmm}},\code{\link[coda]{mcmc}},
  \code{\link[coda]{mcmc.list}}, \code{\link{mcmc.diagnostics.ergmm}} }
\examples{
library(coda)
data(sampson)
monks.fit<-ergmm(samplike~latent(d=2,G=3))
monks.fit.mcmc<-as.mcmc(monks.fit)
plot(monks.fit.mcmc)
raftery.diag(monks.fit.mcmc)
}
\keyword{graphs}
\keyword{debugging}
\keyword{manip}


