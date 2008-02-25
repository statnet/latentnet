\name{plot.ergmm}
\alias{plot.ergmm}
\title{Plotting Method for class ERGMM}
\usage{
\method{plot}{ergmm}(x, \dots, vertex.cex=1, vertex.sides=16*ceiling(sqrt(vertex.cex)),
                          what="mkl",
                          main = NULL, xlab=NULL, ylab=NULL, xlim=NULL,ylim=NULL,
                          object.scale=formals(plot.network.default)$object.scale,
                          pad=formals(plot.network.default)$pad,
                          cluster.col=c("red","green","blue","cyan","magenta","orange","yellow","purple"),
                          vertex.col=NULL, print.formula=TRUE,
                          edge.col=8,
                          pie = FALSE,
                          labels=FALSE,
                          rand.eff=NULL,
                          plot.means=TRUE,plot.vars=TRUE,suppress.axes=FALSE,
                          jitter1D=1,curve1D=TRUE,suppress.center=FALSE,density.par=list())

}
\description{
  \code{\link{plot.ergmm}} is the plotting method for \code{\link[=ergmm.object]{ergmm}} objects.  For
  latent models, this plots the minimum Kullback-Leibler positions by 
  default.  The maximum likelihood, posterior mean, posterior mode, or a
  particular iteration's or configuration's positions
  can be used instead, or pie charts of the posterior
  probabilities of cluster membership can be shown. See
  \code{\link{ergmm}} for more information on how to fit these models.

  At this time, no plotting non-latent-space model fits is not supported.
}
\arguments{
  \item{x}{an R object of class \code{\link[=ergmm.object]{ergmm}}.
    See documentation for \code{\link{ergmm}}.}

  \item{what}{Character vector, integer, or an object of class
    \code{\link[=ergmm.par.object]{ergmm.par}} that specifies
    the point estimates to be used. Can be one of the follwoing:
    \describe{
      \item{\code{"mkl"}}{This is the defult. Plots the Minimum Kulblack-Leibler
	divergence values.}
      \item{\code{"start"},\code{"burnin.start"}}{Plots the starting
	configuration.}
      \item{\code{"sampling.start"}}{Plots the starting configuration of
	the sampling phase (the last burnin configuration).}
      \item{\code{"mle"}}{Plots the maximum likelihood estimates. Random effects
	are treated as fixed.}
      \item{\code{"pmean"}}{Plots the posterior means.}
      \item{\code{"pmode"}}{Plots the conditional posterior mode.}
      \item{\code{"cloud"}}{Plots the ``cloud'' of latent space position draws,
    with their cluster colors.}
      \item{\code{"density"}}{Plots density and contours of the posterior
	latent positions, and, in cluster models, each cluster.}
      \item{\code{\link[=ergmm.par.object]{ergmm.par}}}{Plots the configuration contained in \code{ergmm.fit}.}
      \item{integer}{Plots the configuration of \code{ergm.fit}'th
	stored MCMC iteration.}
    }
  }
  \item{pie}{For latent clustering models, each node is drawn as a pie
    chart representing the probabilities of cluster membership.}
  \item{rand.eff}{A character vector selecting "sender", "receiver",
    "sociality", or "total" random effects. Each vertex is scaled such that its area is proportional
    to the odds ratio due to its selected random effect. }
  \item{plot.means}{Whether cluster means are plotted for latent cluster
    models. The "+" character is used. Defaults to \code{TRUE}.}
  \item{plot.vars}{Whether circles with radius equal to the square root
    of posterior latent or intracluster variance estimates are
    plotted. Defaults to \code{TRUE}.}
  \item{suppress.axes}{Whether axes should \emph{not} be drawn. Defaults
    to \code{FALSE}. (Axes are drawn.)}
  \item{jitter1D}{For 1D latent space fits, it often helps to jitter the
    positions for visualization. This option controls the amount of
    jitter.}
  \item{curve1D}{Controls whether the edges in 1D latent space fits are
    plotted as curves. Defaults to \code{TRUE}.}
  \item{suppress.center}{Suppresses the plotting of "+" at the
    origin. Defaults to \code{FALSE}.}
  \item{cluster.col}{A vector of colors used to distinguish clusters in
    a latent cluster model.}
  \item{main,vertex.cex, vertex.col,
    xlim,ylim,vertex.sides,object.scale,pad,
    edge.col,xlab,ylab}{Arguments passed to
    \code{\link[pkg:network]{plot.network}}, whose defaults differ from
    those of \code{\link[pkg:network]{plot.network}}.}
  \item{labels}{Whether vertex labels should be displayed. Defaults to
    \code{FALSE}.}
  \item{print.formula}{Whether the formula based on which the \code{x}
    was fitted should be printed under the main title. Defaults to \code{TRUE}.}
  \item{density.par}{A list of optional parameters for density plots:
    \describe{
      \item{\code{totaldens}}{Whether the overal density of latent space
	positions should be plotted. Defaults to \code{TRUE}.}
      \item{\code{subdens}}{Whether the densities of latent space
	positions broken down by cluster should be plotted. Defaults to
	\code{TRUE}.}
      \item{\code{mfrow}}{When plotting multiple clusters' densities,
	passed to \code{\link{par}}}
    }
  }
  
  \item{\dots}{Other optional arguments passed to the \code{\link[pkg:network]{plot.network}} function.}
}

\details{
  Plots the results of an ergmm fit.
  
  More information can be found by looking at the documentation of
  \code{\link{ergmm}}.
}
\value{NULL

}
\seealso{\code{\link{ergmm}},\code{\link[=ergmm.object]{ergmm (object)}}, \code{network}, \code{\link[network]{plot.network}}, \code{\link{plot}}}
\examples{
\dontrun{
#
# Using Sampson's Monk data, lets fit a 
# simple latent position model
#
data(sampson)
#
# Using Sampson's Monk data, lets fit a latent clustering random effects model
#
samp.fit <- ergmm(samplike ~ latent(d=2, G=3)+rreceiver)
#
# See if we have convergence in the MCMC
mcmc.diagnostics(samp.fit)
#
# Plot the resulting fit.
#
plot(samp.fit,labels=TRUE,rand.eff="receiver")
plot(samp.fit,pie=TRUE,rand.eff="receiver")
plot(samp.fit,what="pmean",rand.eff="receiver")
plot(samp.fit,what="cloud",rand.eff="receiver")
plot(samp.fit,what="density",rand.eff="receiver")
}
}
\keyword{graphs}
\keyword{hplot}