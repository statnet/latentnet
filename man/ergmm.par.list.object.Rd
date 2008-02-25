\name{ergmm.par.list.object}
\alias{ergmm.par.list.object}
\alias{ergmm.par.list}
\alias{as.ergmm.par.list}
\alias{as.mcmc.list.ergmm.par.list}
\alias{[.ergmm.par.list}
\alias{[[.ergmm.par.list}
\alias{$.ergmm.par.list}
\alias{length.ergmm.par.list}
\alias{stack.ergmm.par.list.list}
\alias{unstack.ergmm.par.list}
\alias{del.iteration}
\title{A List of ERGMM Parameter Configuration}
\description{
  A class \code{\link[=ergmm.par.list.object]{ergmm.par.list}} to represent a
  series of parameter configurations for the same exponential random graph
  mixed model.
}
\details{
  An \code{ergmm.par} object is essentially a named list of lists of parameter
  values. Individual parameters can be accessed with \code{[[} (with
  a character index) or \code{$} operators. \code{[[} operator with a
  numeric or integer index returns an
  \code{\link[=ergmm.par.object]{ergmm.par}} object with the the
  configuration with that index. \code{[} operator given a numeric
  vector returns a \code{ergmm.par.list} object with the subset of
  configurations with the indices given.
  
  The structure of \code{ergmm.par} is derived from the structure of
  \code{\link[=ergmm.par.object]{ergmm.par}}, with each entry having an
  additional dimension (always the first one), indexed by
  configuration. That is, scalars become vectors, vectors become
  matrixes with the original vectors in rows, and matrixes become
  3-dimensional arrays, with the original matrixes indexed by their
  first dimension.
}
\seealso{
  \code{\link{ergmm.par}}, \code{\link{ergmm}}
}
\author{Pavel N. Krivitsky \email{pavel@stat.washington.edu}}
\keyword{graphs}
\keyword{utilities}
\keyword{methods}
\keyword{list}
\keyword{manip}