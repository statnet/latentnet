%  File man/ergmm.prior.Rd in package latentnet, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergmm.prior}
\alias{ergmm.prior}
\title{ Auxiliary for Setting the ERGMM Prior }
\description{
  Auxiliary function as user interface for \code{\link{ergmm}} prior
  specification. Typically only used when calling \code{\link{ergmm}}. It is
  used to supply the parameters of the prior distribution of the model,
  to overwrite those specified in the model formula, and to supply
  miscellaneous prior parameters.
}
\usage{
ergmm.prior(..., adjust.beta.var = TRUE)
}
\arguments{
  \item{\dots}{ Prior distribution parameters. See \link{terms.ergmm}
    for more information. }
  \item{adjust.beta.var}{
    A shortcut: whether the prior variance for each covariate
    coefficient should be divided by the mean square of that
    covariate. This adjustment affects those variances specified in the
    formula or by default, but not those specified through the
    \code{prior=} argument.
  }
}

\value{
  A list with the
  arguments as elements.
}

\seealso{\code{\link{ergmm}}, \code{\link{terms.ergmm}}}

\keyword{ graphs }
\keyword{ models }
