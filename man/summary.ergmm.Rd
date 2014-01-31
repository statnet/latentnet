\name{summary.ergmm}
\alias{summary.ergmm}
\alias{print.summary.ergmm}
\alias{summary.ergmm.object}
\alias{bic.ergmm}
\title{ ERGMM Fit Summaries }
\description{
  \code{summary.ergmm} prodcues a summary of an \code{\link[=ergmm.object]{ergmm}} object, including
  point estimates, standard errors, and BIC calculation.
}
\usage{
\method{summary}{ergmm}(object, point.est=c(
                                     if(!is.null(object[["mle"]])) "mle",
                                     if(!is.null(object[["sample"]])) c("pmean","mkl")
                                     ), quantiles=c(.025,.975),se="mle"\%in\%point.est,\dots)

bic.ergmm(object)
}

\arguments{
  \item{object}{ An \code{\link[=ergmm.object]{ergmm}} object to be summarized.}
  \item{point.est}{ Point estimates to compute: a character vector with
    some subset of \code{"mle"}, \code{"pmean"},
  \code{"mkl"}, and \code{"pmode"}. Defaults
  to a concatenation of \code{"mle"} (if fit), \code{"pmean"}, and
  \code{"mkl"} (if MCMC was run).}
  \item{quantiles}{ Posterior quantiles (credible intervals) to compute. }
  \item{se}{ Whether to compute standard errors. Defaults to \code{TRUE}
  if MLE was fit.}
  \item{\dots}{Additional arguments. Currently unused.}
}
\value{
  For \code{summary}, an object of class
  \code{\link[=summary.ergmm.object]{summary.ergmm}}. A print method is
  available.

  The BICs are available as the element "bic" of the object
  returned. Note that BIC computed for the random effects models may not
  be correct.

  \code{bic.ergmm} returns the BIC for the model directly.
}
\seealso{ \code{\link{ergmm.object}},
  \code{\link{ergmm}}
}

\references{ Chris Fraley and Adrian E. Raftery
  (2002). \emph{Model-based clustering, discriminant analysis, and
  density estimation}. Journal of the American Statistical Association,
  97(458), 611-631.
}

\examples{
\donttest{
data(sampson)
# Fit the model for cluster sizes 1 through 4:
fits<-list(
           ergmm(samplike~euclidean(d=2,G=1)),
           ergmm(samplike~euclidean(d=2,G=2)),
           ergmm(samplike~euclidean(d=2,G=3)),
           ergmm(samplike~euclidean(d=2,G=4))
           )

\dontrun{
# Optionally, plot all fits.
lapply(fits,plot)
}

# Compute the BICs for the fits and plot them:
(bics<-reshape(
    as.data.frame(t(sapply(fits,
                           function(x)c(G=x$model$G,unlist(bic.ergmm(x))[c("Y","Z","overall")])))),
    list(c("Y","Z","overall")),idvar="G",v.names="BIC",timevar="Component",
    times=c("likelihood","clustering","overall"),direction="long"
    ))

with(bics,interaction.plot(G,Component,BIC,type="b",xlab="Clusters", ylab="BIC"))

# Summarize and plot whichever fit has the lowest overall BIC:
bestG<-with(bics[bics$Component=="overall",],G[which.min(BIC)])
summary(fits[[bestG]])
plot(fits[[bestG]])
}
}

\keyword{graphs}
\keyword{models}
\keyword{print}