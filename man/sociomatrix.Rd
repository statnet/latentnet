\name{sociomatrix}
\alias{sociomatrix}
\title{Return the sociomatrix corresponding to a network object}
\usage{
sociomatrix(object, matrix.type="adjacency", attrname = NULL)
}
\description{
  Return the sociomatrix corresponding to a network object as a matrix.
  See \code{\link{ergmm}} for more information.
}
\arguments{
  \item{object}{a \code{network} object. See documentation for the \code{network}
   package.}
  \item{matrix.type}{one of \code{"adjacency"}, \code{"edgelist"}, \code{"incidence"}}
  \item{attrname}{an optional network attribute to return as an edge weight instead of the binary value.}
}

\value{\code{sociomatrix} returns a matrix of zeros and one
(or the attribute value) indicating the presence of a tie from i to j.
}
\seealso{ergmm, network}
\examples{

#
data(sampson)
#
# test the sociomatrix function
#
sociomatrix(samplike)
}
\keyword{graphs}
\keyword{internal}