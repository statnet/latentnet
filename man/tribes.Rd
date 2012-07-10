\name{tribes}
\alias{tribes}
\docType{data}
\title{Read Highland Tribes}
\description{
  A network of political alliances and enmities
  among the 16 Gahuku-Gama sub-tribes of Eastern Central Highlands of New Guinea,
  documented by Read (1954).
}
\usage{data(tribes)}
\format{
  An undirected \code{\link[network]{network}} object with no loops, having the following attributes:
  \describe{
    \item{\code{\%v\% "vertex.names"}}{Character attribute with names of tribes.}
    \item{\code{\%e\% "pos"}}{Logical attribute indicating an
      alliance relationship.}
    \item{\code{\%e\% "neg"}}{Logical attribute indicating a hostile
      relationship ("rova").}
    \item{\code{\%e\% "sign"}}{Numeric attribute coding -1 for enmity, 0
      for no relationship, and 1 for alliance.}
    \item{\code{\%e\% "sign.012"}}{Numeric attribute coding 0 for enmity, 1
      for no relationship, and 2 for alliance.}
  }
  Because of limitations of \code{\link[network]{network}} objects, the object
  itself is a complete graph, and is thus meaningless if used directly
  or plotted.
}
\details{
  This network shows 3 clusters.
}
\source{
  \url{http://vlado.fmf.uni-lj.si/pub/networks/data/UciNet/UciData.htm#gama},
  with corrections from Read (1954).
}
\references{
  Taken from UCINET IV, which cites the following:
  Hage P. and Harary F. (1983). Structural models in anthropology. Cambridge: Cambridge University Press. (See p 56-60).
  Read K. (1954). Cultures of the central highlands, New Guinea. Southwestern Journal of Anthropology, 10, 1-43. 
}
\examples{
\donttest{
data(tribes)
# Only model positive ties:
tribes.fit<-ergmm(tribes~euclidean(d=2,G=3),response="sign.012",family="binomial",fam.par=list(trials=2))
# Edge color must be set manually, for green ties to represent alliance
# and for red ties to represent enmity.
plot(tribes.fit,edge.col=as.matrix(tribes,"pos",m="a")*3+as.matrix(tribes,"neg",m="a")*2,pie=TRUE)
# Model both positive and negative ties:
tribes.fit3<-ergmm(tribes~euclidean(d=2,G=3),response="sign.012",family="binomial.logit",fam.par=list(trials=2))
# Edge color must be set manually, for green ties to represent alliance
# and for red ties to represent enmity.
plot(tribes.fit3,edge.col=as.matrix(tribes,"pos",m="a")*3+as.matrix(tribes,"neg",m="a")*2,pie=TRUE)
}
}
\keyword{datasets}
\keyword{multivariate}
\keyword{cluster}
\keyword{graphs}