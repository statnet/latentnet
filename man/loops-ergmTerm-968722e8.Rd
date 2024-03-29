% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmm.fixed.R
\name{loops-ergmTerm}
\alias{loops-ergmTerm}
\alias{InitErgmTerm.loops}
\title{Self-loops}
\usage{
# binary: loops(mean=0, var=9)

# valued: loops(mean=0, var=9)
}
\arguments{
\item{mean, var}{prior mean and variance.}
}
\description{
Effect of the dyad being a self-loop (i.e., \eqn{(i,i)}).

\strong{Important:} This term works in \CRANpkg{latentnet}'s \code{\link[=ergmm]{ergmm()}} only. Using it in \code{\link[=ergm]{ergm()}} will result in an error.
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "loops", "subsection")}
}
\concept{directed}
\concept{dyad-independent}
\concept{latent}
\concept{undirected}
