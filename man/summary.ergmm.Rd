\name{summary.ergmm}
\alias{summary.ergmm}
\alias{print.summary.ergmm}
\alias{summary.ergmm.object}
\alias{bic.ergmm}
\title{ ERGMM Fit Summaries }
\description{
  \code{summary.ergmm} prodcues a summary of an \code{\link[=ergmm.object]{ergmm}} object, including
  point estimates, standard errors, and BIC calculation.
}
\usage{
\method{summary}{ergmm}(object, point.est = c("pmean", "mkl"),
                        quantiles = c(0.025, 0.975), se = FALSE, \dots)
}

\arguments{
  \item{object}{ An \code{\link[=ergmm.object]{ergmm}} object to be summarized.}
  \item{point.est}{ Point estimates to compute: a character vector with
    some subset of the following: "mle", "pmean","mkl","pmode". }
  \item{quantiles}{ Posterior quantiles (credible intervals) to compute. }
  \item{se}{ Whether to compute standard errors. }
  \item{\dots}{Additional arguments. Currently unused.}
}
\value{
  An object of class
  \code{\link[=summary.ergmm.object]{summary.ergmm}}. A print method is
  available.

  The BICs are available as the element "bic" of the object
  returned. Note that BIC computed for the random effects models may not be correct.
}
\seealso{ \code{\link{ergmm.object}},
  \code{\link{ergmm}}
}
\examples{
data(sampson)
samp.fit <- ergmm(samplike ~ euclidean(d=2, G=3)+rreceiver)
summary(samp.fit)
}
\keyword{graphs}
\keyword{models}
\keyword{print}