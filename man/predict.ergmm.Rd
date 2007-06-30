\name{predict.ergmm}
\alias{predict.ergmm}
\title{ Predicted Dyad Values for an ERGMM. }
\description{
  Returns a matrix of expected dyad values based on an ERGMM fit.
}
\usage{
\method{predict}{ergmm}(ergmm.fit, which.par = "post")
}
\arguments{
  \item{ergmm.fit}{ An obejct of class \code{\link[=ergmm.object]{ergmm}}. }
  \item{which.par}{ One of "mkl", "start", "mle", "pmean", "mkl",
    "pmode", "post", an index of the iteration to use, or a
    \code{\link{ergmm.par}} object, for the set of parameters based on
    which the prediction is made. An exception is "post", which computes
    the expected dyad values integrated over the posterior. }
}
\value{
  A sociomatrix of predicted values.
}
\author{ Pavel N. Krivitsky \email{pavel@stat.washington.edu}}
\seealso{ \code{\link[=ergmm.object]{ergmm}} }
\examples{
data(sampson)
monks.fit<-ergmm(samplike~latent(d=2,G=3))
heatmap(predict(monks.fit))
}
\keyword{ graphs}

