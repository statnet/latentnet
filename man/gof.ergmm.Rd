\name{gof}
\alias{gof.ergmm}
\alias{gof.default}
\alias{gof.formula}
\alias{gof.ergmm}
\alias{gof}
\title{Conduct Goodness-of-Fit Diagnostics on a Exponential Family
  Random Graph Mixed Model Fit}
\description{
  \code{\link{gof}} calculates \eqn{p}-values for geodesic
  distance, degree, and reachability summaries to
  diagnose the goodness-of-fit of exponential family random graph
  mixed models.  See \code{\link{ergmm}} for more information on these models.
}
\usage{
\method{gof}{ergmm}(object, \dots,
         nsim=100,
         GOF=~idegree+odegree+distance,
         verbose=FALSE)
}
\arguments{
  \item{object}{an \code{\link[=ergmm.object]{ergmm}} object (returned by
    \code{\link{ergmm}}).}
  \item{nsim}{The number of simulations to use for the MCMC \eqn{p}-values.
    This is the size of the sample of graphs to be randomly drawn
    from the distribution specified by the object on the set of all graphs.
  }
  \item{GOF}{formula; an \R formula object, of the form
    \code{~ <model terms>} specifying the
    statistics to use to diagnosis the goodness-of-fit of the model. 
    They do not need to be in the model formula specified in
    \code{formula}, and typically are not.
    Examples are the degree distribution ("degree"), 
    minimum geodesic distances ("dist"), and shared partner distributions
    ("espartners" and "dspartners").
    For the details on the possible
    \code{<model terms>}, see \code{\link[ergm]{ergm-terms}}.}
  \item{verbose}{Provide verbose information on the progress of the
    simulation.}
  \item{\dots}{Additional arguments, to be passed to lower-level functions
  in the future.}
}

\details{
  A sample of graphs is randomly drawn from the posterior of the
  \code{\link{ergmm}}.

  A plot of the summary measures is plotted.
  More information can be found by looking at the documentation of
  \code{\link[ergm]{ergm}}.
}
\value{\code{\link{gof}} and \code{\link{gof.ergmm}}
  return an object of class \code{gofobject}.
  This is a list of the tables of statistics and \eqn{p}-values.
  This is typically plotted using \code{\link[ergm]{plot.gofobject}}.
}

\seealso{\code{\link{ergmm}},
  \code{\link[=ergmm.object]{ergmm (object)}},
  \code{\link[ergm]{ergm}}, \code{network},
  \code{\link{simulate.ergmm}}, \code{\link[ergm]{plot.gofobject}}}
\examples{
\donttest{
#
data(sampson)
#
# test the gof.ergm function
#
samplike.fit <- ergmm(samplike ~ euclidean(d=2,G=3),
                      control=ergmm.control(burnin=1000,interval=5))
samplike.fit
summary(samplike.fit)

#
# Plot the probabilities first
#
monks.gof <- gof(samplike.fit)
monks.gof
#
# Place all three on the same page
# with nice margins
#
par(mfrow=c(1,3))
par(oma=c(0.5,2,1,0.5))
#
plot(monks.gof)
#
# And now the odds 
#
plot(monks.gof, plotlogodds=TRUE)
}
}
\keyword{models}
