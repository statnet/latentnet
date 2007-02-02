\name{ergmm}
\alias{ergmm}
\alias{latent}
\alias{latentcluster}
\title{Fit a Latent Space Random Graph Model}
\description{
  \code{\link{ergmm}} is used to fit latent space and latent space cluster random network models, 
    as described in Hoff, Raftery and Handcock (2002),
    Handcock, Raftery and Tantrum (2005), and the article we are working on.
  \code{\link{ergmm}} can return either a Bayesian model fit
  or the two-stage MLE.
}
\usage{
ergmm(formula, response=NULL, family="Bernoulli",fam.par=NULL,
      control=ergmm.control(), user.start=NULL, prior=list(),
      skipMCMC=FALSE, skipProcrustes=FALSE, skipMBC=FALSE,
      Z.ref=NULL, Z.K.ref=NULL, store.burnin=FALSE, randseed=NULL,
      verbose=FALSE, orthogonalize=FALSE)
}
\arguments{
  \item{formula}{An \R formula object, of the form
    \code{g ~ <term 1> + <term 2> ...},
    where \code{g} is a network object or a matrix that can be coerced to a
    network object, and \code{<term 1>}, \code{<term 2>}, etc., are each 
    terms for the model. See \code{\link{terms.ergmm}} for the terms
    that can be fitted.
    To create a network object in \R, use the \code{network} function,
    then add nodal attributes to it using \code{set.vertex.attribute} 
    if necessary.}
  \item{response}{An optional edge attribute that serves as the response
    variable. By default, presence (1) or absence (0) of an edge in
    \code{g} is used.}
  \item{family}{A character vector that is one of "Bernoulli" (the default),
    "binomial", or "Poisson", specifying the conditional distribution of
    each edge value. See \code{Details} for more information.}
  \item{fam.par}{For those families that require additional parameters,
    a list.}
  \item{control}{The MCMC parameters that do not affect the posterior
    distribution such as the sample size, the proposal variances, and
    tuning parameters, in the
    form of a named list. See \code{\link{ergmm.control}} for more
    information and defaults.}
  \item{user.start}{An optional initial configuration parameters for
    MCMC in the form of a list. By default, posterior mode conditioned on cluster assignments
    is used. It is permitted to only supply some of the parameters of a
    configuration. If this is done, the remaining paramters are fitted
    conditional on those supplied.}
  \item{prior}{The prior parameters for the model being fitted in the
    form of a named list. See \link{terms.ergmm} for the names to use.
    If given, will override those given in the
    formula terms, making it useful as a convenient way to store and reproduce a
    prior distribution. The list or prior parameters can also be
    extracted from an \link[=ergmm.object]{ERGMM fit object}. Defaults to an empty list (that is, all prior
    parameters are set in the formula terms).}
  \item{skipMCMC}{Skip the MCMC sampling and return only the approximate
    posterior mode. Defaults to \code{FALSE}.}
  \item{skipProcrustes}{Skip the Procrustes and label-switching of the
    MCMC sample. Defaults to \code{FALSE}.}
  \item{skipMBC}{Skip the Bayesian model-based clustering run to find
    posterior cluster means and variances consistent with MKL positional estimates. Defaults to \code{FALSE}.}
  \item{Z.ref}{If given, used as a reference for Procrustes analysis.}
  \item{Z.K.ref}{If given, used as a reference for label-switching.}
  \item{store.burnin}{If \code{TRUE}, the MCMC samples from the burnin
    period are stored (though not used in the inference). This option is
    mainly useful for diagnosing slow convergence. Defaults to
    \code{FALSE}.}
 \item{randseed}{If supplied, random number integer seed.}
 \item{verbose}{If this is \code{TRUE}, causes information to be
   printed out about the progress of the fitting, particularly initial
   value generation. Higher values lead to more verbosity.}
 \item{orthogonalize}{\emph{EXPERIMENTAL.} Performs Gramm-Schmidt
   orthogonalization on covariates passed to the model to reduce
   dependence. This capacity may be obsoleted in the future by the use
   of correlated proposals.}
}
\value{
  \code{\link{ergmm}} returns an object of class \code{\link[=ergmm.object]{ergmm}}
  containing the information about the posterior.
}
\seealso{network, set.vertex.attributes, set.network.attributes,
         summary.ergmm, print.ergmm,\link{terms.ergmm}}
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
#
# See http://www.csde.washington.edu/statnet
# for examples
#
# load the Florentine marriage data matrix
#
data(flo)
#
# attach the sociomatrix for the Florentine marriage data
# This is not yet a network object.
#
flo
#
# Create a network object out of the adjacency matrix
#
flomarriage <- network(flo,directed=FALSE)
flomarriage
#
# print out the sociomatrix for the Florentine marriage data
#
sociomatrix(flomarriage)
#
# create a vector indicating the wealth of each family (in thousands of lira) 
# and add it as a covariate to the network object
#
flomarriage <- set.vertex.attribute(flomarriage,"wealth",
  c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3))
flomarriage
#
# create a plot of the social network
#
plot(flomarriage)
#
# now make the vertex size proportional to their wealth
#
plot(flomarriage, vertex.cex="wealth", main="Marriage Ties")
#
# Use 'data(package = "latentnet")' to list the data sets in a
#
data(package="latentnet")
#
# Using Sampson's Monk data, lets fit a 
# simple latent position model
#
data(sampson)
#
# Get the group labels
#
group <- get.vertex.attribute(samplike,"group")
samp.labs <- substr(group,1,1)
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
\dontrun{
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
#
# Add contours
#
ergmm.add.contours(samp.fit,nlevels=8,lwd=2)
points(samp.fit$Z.mkl,pch=19,col=samp.fit$class)
#
# Try a covariate on the group
#
samegroup <- outer(group, group, "==")
diag(samegroup) <- 0
samp.fit <- ergmm(samplike ~ latentcov(samegroup) + latent(k=2))
summary(samp.fit)
}
}
\keyword{graphs}
