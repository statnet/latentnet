\name{ergmm.drawpie}
\alias{ergmm.drawpie}
\title{
  Draw a pie chart at a specified location.
}
\description{
  Used by \code{\link{plot.ergmm}} to draw pie charts to visualize soft
  clusterings when \code{pie=TRUE}. Exported as a courtesy to dependent
  packages.
}
\usage{
ergmm.drawpie(center, radius, probs, n = 50, cols = 1:length(probs), ...)
}
\arguments{
  \item{center}{
    A numeric vector of length 2, specifying the horizontal and the
    vertical coordinates of its center.
}
  \item{radius}{
    Radius of the pie chart.
}
  \item{probs}{
    A vector of probabilities/weights of each sector; they do not have
    to sum to 1.
}
  \item{n}{
    Number of points to use to approximate the "circle".
}
  \item{cols}{
    A vector of colors to use for the sectors.
}
  \item{\dots}{
    Additional arguments, currently unused.
}
}
\author{
  See COPYRIGHT.
}

\seealso{
  plot.ergmm
}
\examples{
plot(c(0,sum(1:11))*2,c(-10,10),type="n",asp=1)
for(i in 1:10) ergmm.drawpie(c(sum(1:i)*2,0), radius=i, probs=1:(i+1))
}
\keyword{ graphs }
