\name{terms.ergmm}
\alias{terms.ergmm}
\alias{ergmm.terms}
\title{Model Terms for Latent Space Random Graph Model}
\description{
  Model terms that can be used in an \code{\link{ergmm}} formula.
}
\section{Model Terms}{
  The \code{\link{latentnet}} package itself allows only
  dyad-independent terms. In the formula for the model, the model terms are various function-like
  calls, some of which require arguments, separated by \code{+} signs.
  
  \emph{Random Effects}\cr
  \describe{
    \item{\code{latent(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
  mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)}}{\emph{Latent position model term, with
	optional clustering.}
      Adds a term to the model equal to the negative Eucledean distance
      \eqn{-||Z_i-Z_j||}{-dist(Z[i],Z[j])}, where \eqn{Z_i}{Z[i]} and \eqn{Z_j}{Z[j]}
      are the positions of their respective actors in an unobserved social
      space. These positions may optionally have a finite spherical
      Gaussian mixture clustering structure.
      The parameters are as follows:
      \describe{
	\item{\code{d}}{The dimension of the latent space.}
	\item{\code{G}}{The number of groups (0 for no clustering).}
	\item{\code{var.mul}}{In the absence of \code{var}, this
	  argument will be used as a scaling factor for a
	  function of average cluster size and latent space dimension to
	  set \code{var}. To set it in the \code{prior} argument to \code{\link{ergmm}}, use \code{Z.var.mul}.}
	\item{\code{var}}{If given, the scale
	  parameter for the scale-inverse-chi-squared prior distribution of the within-cluster
	  variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
	  \code{Z.var}.}
	\item{\code{var.df.mul}}{In the absence of \code{var.df}, this
	  argument is the multiplier for the square root of average
	  cluster size, which serves in place of \code{var.df}. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
	  \code{Z.var.df.mul}.}
	\item{\code{var.df}}{The degrees of freedom
	  parameter for the scale-inverse-chi-squared prior distribution of the within-cluster
	  variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
	  \code{Z.var.df}.}
	\item{\code{mean.var.mul}}{In the absence of \code{mean.var},
	  the multiplier for a function of number of vertices and latent space
	  dimension to set \code{mean.var}. To set it in the
	  \code{prior} argument to \code{\link{ergmm}}, use \code{Z.mean.var.mul}.}
	\item{\code{mean.var}}{The variance of
	  the spherical Gaussian prior distribution of the cluster means. To set it in the
	  \code{prior} argument to \code{\link{ergmm}}, use
	  \code{Z.mean.var}.}
	\item{\code{pK.mul}}{In the absence of \code{pK}, this argument
	  is the multiplier for the square root of the average cluster size,
	  which is used as \code{pK}. To set it in
	  the \code{prior} argument to \code{\link{ergmm}}, use \code{Z.pK}.}
	\item{\code{pK}}{The parameter of the Dirichilet prior
	  distribution of cluster assignment probabilities. To set it in
	  the \code{prior} argument to \code{\link{ergmm}}, use \code{Z.pK}.}
      }
    }
  }
  
  \emph{Fixed Effects}\cr
  Each coefficient for a fixed effect covariate has a normal prior whose
  mean and variance are set by the \code{mean} and \code{var} parameters
  of the term. For those formula terms that add more than one covariate,
  a vector can be given for mean and variance. If not, the vectors given
  will be repeated until the needed length is reached.
  \describe{
    \item{\code{latentcov(x, attrname=NULL, mean=0, var=9)}}{\emph{Edge covariates for the
	latent model.}
      \code{x} is either a matrix of
      covariates on each pair of vertices, a network, or an edge attribute on \code{g}; 
      if the latter, optional argument
      \code{attrname} provides the name of the edge attribute to
      use for edge values. \code{latentcov} can be called more
      than once, to model the effects of multiple covariates. Note that
      some covariates can be more conveniently specified using the
      following terms. 
    }
    \item{\code{absdiff(attrname, mean=0, var=9)}}{\emph{Absolute Difference.}
      \code{attrname} is a character string giving the
      name of an attribute in the network's vertex attribute list. This
      term adds a covariate with the value
      \code{abs(attrname(i)-attrname(j))} for all edges.   
    }
    \item{\code{nodematch(attrname, diff=FALSE, mean=0, var=9)}}{\emph{Uniform Homophily
	and Differential Homophily.}  \code{attrname} is a
      character string giving the name of an attribute in the
      network's vertex attribute list.  When \code{diff=FALSE},
      this term adds one covariate with the value
      \code{attrname(i)==attrname(j)}.  When \code{diff=TRUE},
      \eqn{p} covariates are added to the model, where \eqn{p} is the
      number of unique values of the \code{attrname} attribute.
      The \eqn{k}th such covariate has the value \code{attrname(i) == attrname(j) == value(k)}, where
      \code{value(k)} is the \eqn{k}th smallest unique value of the
      \code{attrname} attribute.
    }
    \item{\code{sendercov(attrname, force.factor=FALSE, mean=0, var=9)}}{\emph{Sender covariate effect.}
      \code{attrname} is a character string giving the name of an
      attribute in the network's vertex attribute list.
      If the attribute is numeric, This term adds one covariate
      to the model equaling \code{attrname(i)}. If the attribute is not
      numeric or \code{force.factor==TRUE}, this term adds \eqn{p-1}
      covariates to the model,
      where \eqn{p} is the number of unique values of \code{attrname}.
      The \eqn{k}th such covariate has the value \code{attrname(i) == value(k+1)}, where
      \code{value(k)} is the \eqn{k}th smallest unique value of the
      \code{attrname} attribute. This term only makes
      sense if \code{g} is directed.}
    \item{\code{receivercov(attrname, force.factor=FALSE, mean=0, var=9)}}{\emph{Receiver covariate effect.}
      \code{attrname} is a character string giving the name of an
      attribute in the network's vertex attribute list.
      If the attribute is numeric, This term adds one covariate
      to the model equaling \code{attrname(j)}. If the attribute is not
      numeric or \code{force.factor==TRUE}, this term adds \eqn{p-1}
      covariates to the model,
      where \eqn{p} is the number of unique values of \code{attrname}.
      The \eqn{k}th such covariate has the value \code{attrname(j) == value(k+1)}, where
      \code{value(k)} is the \eqn{k}th smallest unique value of the
      \code{attrname} attribute. This term only makes
      sense if \code{g} is directed.}
    \item{\code{socialitycov(attrname, force.factor=FALSE, mean=0, var=9)}}{\emph{Sociality covariate effect.}
      \code{attrname} is a character string giving the name of an
      attribute in the network's vertex attribute list.
      If the attribute is numeric, This term adds one covariate
      to the model equaling \code{attrname(i)+attrname(j)}. If the attribute is not
      numeric or \code{force.factor==TRUE}, this term adds \eqn{p-1}
      covariates to the model,
      where \eqn{p} is the number of unique values of \code{attrname}.
      The \eqn{k}th such covariate has the value \code{attrname(i) ==
	value(k+1) + attrname(j) == value(k+1)}, where
      \code{value(k)} is the \eqn{k}th smallest unique value of the
      \code{attrname} attribute. This term makes sense whether or not
      \code{g} is directed.}
  }
}
\seealso{
  \code{\link{ergmm}}
}
\keyword{graphs}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
\keyword{nonparametric}
\keyword{cluster}
