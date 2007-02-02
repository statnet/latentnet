\name{latentnet}
\alias{latentnet}
\title{Latent position and cluster models for networks}
\description{
  The package \code{statnet} is used to fit linear exponential random network
  models, in which the probability of a given network, \eqn{g},
  on a set of nodes is \eqn{\exp(\theta^{T}S(g))/c(\theta)},
  where \eqn{S(g)} is a vector of network statistics, \eqn{\theta}
  is a parameter vector of the same length and \eqn{c(\theta)}
  is the normalizing constant for the distribution.  \code{\link{ergmm}}
  can return either a maximum pseudo-likelihood estimate or
  an approximate MLE based on a Monte Carlo scheme.

  The package \code{\link{latentnet}} adds latent variable models to the traditional 
  exponential random network models. These latent models are
  fit using the \code{\link{ergmm}} function that is also used
  by \code{statnet}. See \url{http://www.csde.washington.edu/statnet}.

  The \code{\link{ergmm}} specifies models via: \code{g ~ <model terms>}
  where \code{g} is a \code{network} object
  For the list of possible \code{<model terms>}, see 
% the \code{Details} section below.
  the \code{Model Terms} section of the \code{\link{ergmm}} function.
}
\usage{
}
\arguments{
\item{\dots}{Specific to the latent variable model terms.}
}

\details{

    The arguments in the \code{\link{ergmm}} function specific to latent 
    variable models
    are \code{latent.control}. See the help page for \code{\link{ergmm}} for the 
    details.

    The result of a latent variable model fit is an \code{\link{ergmm}} object.
    Hence the \code{\link{summary}}, \code{print}, and \code{plot} functions
    apply to the fits.
    The \code{\link{plot.ergmm}} function has many options specific to latent 
    variable models.
% are \code{mle}, \code{pie}, \code{contour},
% \code{contour.colors}.
    See the help page for \code{\link{plot.ergmm}} for the details.
}
\value{\code{\link{ergmm}} returns an object of class 'ergmm' that is a list.
}
\seealso{latent, latentcluster, sna, network, terms.ergmm}
\references{
 Peter D. Hoff,  Adrian E. Raftery and Mark S. Handcock.
 \emph{Latent space approaches to social network analysis.} 
 Journal of the American Statistical Association, Dec 2002, Vol.97, 
 Iss. 460;  pg. 1090-1098.

 Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum.
 \emph{Model-Based Clustering for Social Networks.} 
 Working Paper Number 46, Center for Statistics and the Social Sciences,
 University of Washington, April 2005.
}
\examples{
%\dontshow{require(statnet)}
#
# Using Sampson's Monk data, lets fit a 
# simple latent position model
#
data(sampson)
#
# Get the group labels
samp.labs <- substr(get.vertex.attribute(samplike,"group"),1,1)
#
samp.fit <- ergmm(samplike ~ latent(k=2), burnin=10000,
                 MCMCsamplesize=2000, interval=30)
#
# See if we have convergence in the MCMC
mcmc.diagnostics(samp.fit)
#
# Plot the fit
#
plot(samp.fit,label=samp.labs, vertex.col="group")
#
# Using Sampson's Monk data, lets fit a latent clustering model
#
samp.fit <- ergmm(samplike ~ latentcluster(k=2, ngroups=3), burnin=10000,
                 MCMCsamplesize=2000, interval=30)
#
# See if we have convergence in the MCMC
mcmc.diagnostics(samp.fit)
#
# Lets look at the goodness of fit:
#
plot(samp.fit,label=samp.labs, vertex.col="group")
plot(samp.fit,pie=TRUE,label=samp.labs)
plot(samp.fit,density=c(2,2))
plot(samp.fit,contours=5,contour.color="red")
plot(samp.fit,density=TRUE,drawarrows=TRUE)
ergmm.add.contours(samp.fit,nlevels=8,lwd=2)
points(samp.fit$Z.mkl,pch=19,col=samp.fit$class)
}
\keyword{graphs}
