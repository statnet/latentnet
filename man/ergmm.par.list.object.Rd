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
  \code{[[} operator with a
  numeric or integer index returns a
  list with the the
  configuration with that index. \code{[} operator given a numeric
  vector returns a \code{ergmm.par.list} object with the subset of
  configurations with the indices given.
  
  The structure of \code{ergmm.par.list} is derived from named lists, with each entry having an
  additional dimension (always the first one), indexed by
  configuration. That is, scalars become vectors, vectors become
  matrixes with the original vectors in rows, and matrixes become
  3-dimensional arrays, with the original matrixes indexed by their
  first dimension.
}
\seealso{
  \code{\link{ergmm}}
}
\keyword{graphs}
\keyword{utilities}
\keyword{methods}
\keyword{list}
\keyword{manip}