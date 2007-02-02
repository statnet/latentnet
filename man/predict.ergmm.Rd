\name{predict.ergmm}
\alias{predict.ergmm}
\title{ Predicted Dyad Values for an ERGMM. }
\description{
  Returns a matrix of expected dyad values based on an ERGMM fit.
}
\usage{
\method{predict}{ergmm}(ergmm.fit, theta = "mkl")
}
\arguments{
  \item{ergmm.fit}{ An obejct of class \code{\link[=ergmm.object]{ergmm}}. }
  \item{theta}{ One of "mkl", "start", "mle", "pmean", "mkl", "pmode",
    an index of the iteration to use, or a \code{\link{ergmm.par}} object. }
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

