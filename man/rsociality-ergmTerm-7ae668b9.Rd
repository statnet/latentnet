% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmm.random.R
\name{rsociality-ergmTerm}
\alias{rsociality-ergmTerm}
\alias{InitErgmTerm.rsociality}
\title{Random sociality effect}
\usage{
# binary: rsociality(var=1, var.df=3)

# valued: rsociality(var=1, var.df=3)
}
\arguments{
\item{var}{The scale parameter for the scale-inverse-chi-squared
prior distribution of the sociality effect variance. To set
it in the \code{prior} argument to \code{\link{ergmm}}, use
\code{sociality.var}.}

\item{var.df}{The degrees of freedom parameter for the
scale-inverse-chi-squared prior distribution of the sociality effect
variance. To set it in the \code{prior} argument to
\code{\link{ergmm}}, use \code{sociality.var.df}.}
}
\description{
Adds a random sociality effect to the model, with normal
prior centered around \eqn{0}{0} and a variance that is
estimated.

\strong{Important:} This term works in \CRANpkg{latentnet}'s \code{\link[=ergmm]{ergmm()}} only. Using it in \code{\link[=ergm]{ergm()}} will result in an error.
}
\details{
The following parameters are associated with this term:
\describe{

\item{\code{sociality}}{ Numeric vector of values of each
vertex's random sociality effect.}

\item{\code{sociality.var}}{ Random sociality effect's variance.}

}
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "rsociality", "subsection")}
}
\concept{directed}
\concept{dyad-independent}
\concept{latent}
\concept{undirected}
