% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmm.random.R
\name{rreceiver-ergmTerm}
\alias{rreceiver-ergmTerm}
\alias{InitErgmTerm.rreceiver}
\title{Random receiver effect}
\usage{
# binary: rreceiver(var=1, var.df=3)

# valued: rreceiver(var=1, var.df=3)
}
\arguments{
\item{var}{The scale parameter for the scale-inverse-chi-squared
prior distribution of the receiver effect variance. To set
it in the \code{prior} argument to \code{\link{ergmm}}, use
\code{receiver.var}.}

\item{var.df}{The degrees of freedom parameter for the
scale-inverse-chi-squared prior distribution of the receiver effect
variance. To set it in the \code{prior} argument to
\code{\link{ergmm}}, use \code{receiver.var.df}.}
}
\description{
Adds a random receiver effect to the model, with normal
prior centered around \eqn{0}{0} and a variance that is
estimated. Can only be used on directed networks.

\strong{Important:} This term works in \CRANpkg{latentnet}'s \code{\link[=ergmm]{ergmm()}} only. Using it in \code{\link[=ergm]{ergm()}} will result in an error.
}
\details{
The following parameters are associated with this term:
\describe{

\item{\code{receiver}}{ Numeric vector of values of each
vertex's random receiver effect.}

\item{\code{receiver.var}}{ Random receiver effect's variance.}

}
}
\note{
This term can only be used with directed networks.
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "rreceiver", "subsection")}
}
\concept{directed}
\concept{dyad-independent}
\concept{latent}
