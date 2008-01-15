\name{families.ergmm}
\alias{families.ergmm}
\alias{ergmm.families}
\alias{dlpY.deta.Bernoulli.logit}
\alias{dlpY.deta.binomial.logit}
\alias{dlpY.deta.fs}
\alias{dlpY.deta.Poisson.log}
\alias{lpYc.Bernoulli.logit}
\alias{lpYc.binomial.logit}
\alias{lpYc.fs}
\alias{lpYc.Poisson.log}
\alias{lpY.Bernoulli.logit}
\alias{lpY.binomial.logit}
\alias{lpY.Poisson.log}
\alias{lpY.fs}
\alias{EY.Bernoulli.logit}
\alias{EY.binomial.logit}
\alias{EY.fs}
\alias{EY.Poisson.log}
\alias{mk.lp.optim.fs}
\alias{pY.fs}
\alias{pY.Poisson.log}
\alias{pY.Bernoulli.logit}
\alias{pY.binomial.logit}
\alias{rsm.fs}
\alias{rsm.Poisson.log}
\alias{rsm.binomial.logit}
\alias{rsm.Bernoulli.logit}
\alias{family.IDs}
\alias{family.names}
\alias{fam.par.check}
\title{Edge Weight Distribution Families}
\description{
  Family-link combinations supported by \code{\link{ergmm}}.
}
\section{Family-link combinations}{
  \tabular{rlllll}{
    ID \tab \code{C} name          \tab \R name   \tab Type       \tab Family    \tab Link  \cr
     1 \tab Bernoulli\_logit       \tab Bernoulli.logit \tab Discrete   \tab Bernoulli \tab logit \cr
     2 \tab binomial\_logit        \tab binomial.logit  \tab Discrete   \tab binomial  \tab logit \cr
     3 \tab Poisson\_log           \tab Poisson.log     \tab Discrete   \tab Possion   \tab log   \cr
     4 \tab Bernoulli\_cont\_logit \tab NA              \tab Continuous \tab Bernoulli \tab logit \cr
     5 \tab binomial\_cont\_logit  \tab NA              \tab Continuous \tab binomial  \tab logit \cr
     6 \tab Poisson\_cont\_log     \tab NA              \tab Continuous \tab Possion   \tab log  
   }
}
\details{
  Each supported family has a family of \R functions, of the form
  \code{pY.}-, \code{lpY.}-, \code{EY.}-,
  \code{dlpY.deta.}-, \code{lpYc.}-, \code{rsm.}-, followed by the family's \R name,
  for the respective family's \R name, representing the family's
  likelihood, log-likelihood, expectation, derivative of log-likelihood,
  log-normalizing-constant, and random sociomatrix generation functions.

  On the \code{C} side, similar functions exist, but becuase of static
  typing, are also provided for ``continuous'' versions of those families.
}
\keyword{graphs}
\keyword{models}
\keyword{regression}