\name{summary.ergmm}
\alias{summary.ergmm}
\title{ ERGMM Fit Summaries }
\description{
  \code{summary.ergmm} prodcues a summary of an \code{\link[=ergmm.object]{ergmm}} object, including
  point estimates and standard errors.
}
\usage{
\method{summary}{ergmm}(ergmm.fit, point.est = c("pmean", "mkl"), quantiles = c(0.025, 0.975), se = TRUE)
}

\arguments{
  \item{ergmm.fit}{ An \code{\link[=ergmm.object]{ergmm}} object to be summarized.}
  \item{point.est}{ Point estimates to compute: a character vector with
    some subset of the following: "mle", "pmean","mkl","pmode". }
  \item{quantiles}{ Posterior quantiles (credible intervals) to compute. }
  \item{se}{ Whether to compute standard errors. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  An object of class \code{\link[=summary.ergmm.object]{summary.ergmm}}.
}
\author{ Pavel N. Krivitsky \name{pavel@stat.washington.edu} }
\seealso{ \code{\link{ergmm.object},
\code{\link{ergmm}}
}
\examples{
}
\keyword{ graphs }
