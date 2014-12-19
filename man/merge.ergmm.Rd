%  File man/merge.ergmm.Rd in package latentnet, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{merge.ergmm}
\alias{merge.ergmm}
\title{ Merge two or more replications of ERGMM fits }
\description{
  \code{merge.ergmm} prodcues a \code{\link[=ergmm.object]{ergmm}}
  object containing the combined MCMC output (and derived estimates) of
  several \code{\link[=ergmm.object]{ergmm}} objects produced with the
  same input parameters but different starting values, random seeds, etc..
}
\usage{
\method{merge}{ergmm}(x, y, \dots, verbose = FALSE)
}

\arguments{
  \item{x}{ The first \code{\link[=ergmm.object]{ergmm}} object to be
    merged.}
  \item{y}{ The second \code{\link[=ergmm.object]{ergmm}} object to be merged.}
  \item{\dots}{ Additional \code{\link[=ergmm.object]{ergmm}} objects to
  be merged.}
  \item{verbose}{If \code{TRUE}, marks the progress of merging.}
}
\value{
  An object of class
  \code{\link[=ergmm.object]{ergmm}}. 
}
\seealso{ \code{\link{ergmm.object}},
  \code{\link{ergmm}}
}
\examples{
\donttest{
data(sampson)
# Run two short MCMC-based fits.
samp.fit1 <- ergmm(samplike ~ euclidean(d=2, G=3),
  control=ergmm.control(burnin=1000,interval=10,sample.size=2000))
samp.fit2 <- ergmm(samplike ~ euclidean(d=2, G=3),
  control=ergmm.control(burnin=1000,interval=10,sample.size=2000))

# Combine them, and summarize the result.
samp.fit <-  merge.ergmm(samp.fit1,samp.fit2)
summary(samp.fit)
}
}
\keyword{graphs}
\keyword{models}
