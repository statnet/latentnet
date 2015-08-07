\name{terms.ergmm}
\alias{terms.ergmm}
\alias{ergmm.terms}
\alias{terms-ergmm}
\alias{ergmm-terms}
\alias{dlpY.dZ.bilinear}
\alias{dlpY.dZ.negative.Euclidean}
\alias{dlpY.dZ.negative.Euclidean2}
\alias{dlpY.dZ.fs}
\alias{latent.effect.fs}
\alias{latent.effect.IDs}
\alias{latent.effect.invariances}
\alias{latent.effect.names}
\alias{InitErgmm.intercept}
\alias{InitErgmm.Intercept}
\alias{InitErgmm.1}
\alias{InitErgmm.loops}
\alias{InitErgmm.loopcov}
\alias{InitErgmm.loopfactor}
\alias{InitErgmm.euclidean}
\alias{InitErgmm.bilinear}
\alias{InitErgmm.euclidean2}
\alias{InitErgmm.latentcov}
\alias{InitErgmm.receivercov}
\alias{InitErgmm.rreceiver}
\alias{InitErgmm.rsender}
\alias{InitErgmm.rsociality}
\alias{InitErgmm.sendercov}
\alias{InitErgmm.socialitycov}
\alias{intercept}
\alias{Intercept}
\alias{1}
\alias{loops}
\alias{loopcov}
\alias{loopfactor}
\alias{euclidean}
\alias{bilinear}
\alias{euclidean2}
\alias{latentcov}
\alias{receivercov}
\alias{rreceiver}
\alias{rsender}
\alias{rsociality}
\alias{sendercov}
\alias{socialitycov}
\alias{latent.effect.bilinear}
\alias{latent.effect.negative.Euclidean}
\alias{latent.effect.negative.Euclidean2}

\title{Model Terms for Latent Space Random Graph Model}
\description{
  Model terms that can be used in an \code{\link{ergmm}} formula and
  their parameter names.
}
\section{Model Terms}{
  The \code{\link{latentnet}} package itself allows only
  dyad-independent terms. In the formula for the model, the model terms are various function-like
  calls, some of which require arguments, separated by \code{+} signs.

  \emph{Latent Space Effects}\cr
  \describe{
    \item{\code{euclidean(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
  mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)}}{\emph{(Negative)
  Euclidean distance model term, with
	optional clustering.}
      Adds a term to the model equal to the negative Eucledean distance
      \eqn{-||Z_i-Z_j||}{-dist(Z[i],Z[j])}, where \eqn{Z_i}{Z[i]} and \eqn{Z_j}{Z[j]}
      are the positions of their respective actors in an unobserved social
      space. These positions may optionally have a finite spherical
      Gaussian mixture clustering structure. This term was previously
      called \code{latent} which now fits negative Euclidean latent
      space model with a warning.
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
      The following parameters are associated with this term:
      \describe{
	\item{\code{Z}}{ Numeric matrix with rows being latent space
	  positions.}
	\item{\code{Z.K} (when \eqn{\code{G}>0})}{ Integer vector of cluster assignments. }
	\item{\code{Z.mean} (when \eqn{\code{G}>0})}{ Numeric matrix with rows being cluster means. }
	\item{\code{Z.var}  (when \eqn{\code{G}>0})}{ Depending on the model, either a numeric vector with
	  within-cluster variances or a numeric scalar with the overal latent space variance. }
	\item{\code{Z.pK}  (when \eqn{\code{G}>0})}{ Numeric vector of probabilities of a vertex being
	  in a particular cluster.} 
      }
    }
    \item{\code{bilinear(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
  mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)}}{\emph{Bilinear
  latent model term, with
	optional clustering.}
      Adds a term to the model equal to the inner product of the latent positions:
      \eqn{Z_i \cdot Z_j}{sum(Z[i]*Z[j])}, where \eqn{Z_i}{Z[i]} and \eqn{Z_j}{Z[j]}
      are the positions of their respective actors in an unobserved social
      space. These positions may optionally have a finite spherical
      Gaussian mixture clustering structure. \emph{Note: For a
	bilinear latent space effect, two actors being closer in the
	clustering sense does not necessarily mean that the expected value of
	a tie between them is higher. Thus, a warning is printed when
	this model is combined with clustering.}
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
      The following parameters are associated with this term:
      \describe{
	\item{\code{Z}}{ Numeric matrix with rows being latent space
	  positions.}
	\item{\code{Z.K} (when \eqn{\code{G}>0})}{ Integer vector of cluster assignments. }
	\item{\code{Z.mean} (when \eqn{\code{G}>0})}{ Numeric matrix with rows being cluster means. }
	\item{\code{Z.var}  (when \eqn{\code{G}>0})}{ Depending on the model, either a numeric vector with
	  within-cluster variances or a numeric scalar with the overal latent space variance. }
	\item{\code{Z.pK}  (when \eqn{\code{G}>0})}{ Numeric vector of probabilities of a vertex being
	  in a particular cluster.} 
      }
    }
    \item{\code{euclidean2(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
  mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)}}{\emph{(Negative)
  squared Euclidean distance model term, with
	optional clustering.}
      Adds a term to the model equal to the negative Eucledean distance
      \eqn{-||Z_i-Z_j||^2}{-dist(Z[i],Z[j])^2}, where \eqn{Z_i}{Z[i]} and \eqn{Z_j}{Z[j]}
      are the positions of their respective actors in an unobserved social
      space. These positions may optionally have a finite spherical
      Gaussian mixture clustering structure. This term was previously
      called \code{latent} which now fits negative Euclidean latent
      space model with a warning.
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
      The following parameters are associated with this term:
      \describe{
	\item{\code{Z}}{ Numeric matrix with rows being latent space
	  positions.}
	\item{\code{Z.K} (when \eqn{\code{G}>0})}{ Integer vector of cluster assignments. }
	\item{\code{Z.mean} (when \eqn{\code{G}>0})}{ Numeric matrix with rows being cluster means. }
	\item{\code{Z.var}  (when \eqn{\code{G}>0})}{ Depending on the model, either a numeric vector with
	  within-cluster variances or a numeric scalar with the overal latent space variance. }
	\item{\code{Z.pK}  (when \eqn{\code{G}>0})}{ Numeric vector of probabilities of a vertex being
	  in a particular cluster.} 
      }
    }
  }
  \emph{Actor-specific effects}\cr
  \describe{
    \item{\code{rsender(var=1, var.df=3)}}{\emph{Random sender effect.}
      Adds a random sender effect to the model, with normal prior centered
      around \eqn{0}{0} and a variance that is estimated.
      Can only be used on directed networks.
      The parameters are as follows:
      \describe{
	\item{\code{var}}{The scale
	  parameter for the scale-inverse-chi-squared prior distribution
	  of the sender effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
	  \code{sender.var}.}
	\item{\code{var.df}}{The degrees of freedom
	  parameter for the scale-inverse-chi-squared prior distribution
	  of the sender effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
	  \code{sender.var.df}.
	}
      }
      The following parameters are associated with this term:
      \describe{
	\item{\code{sender}}{ Numeric vector of values of each
	  vertex's random sender effect.}
	\item{\code{sender.var}}{ Random sender effect's variance.}
      }
    }
    \item{\code{rreceiver(var=1, var.df=3)}}{\emph{Random receiver effect.}
      Adds a random receiver effect to the model, with normal prior centered
      around \eqn{0}{0} and a variance that is estimated.
      Can only be used on directed networks.
      The parameters are as follows:
      \describe{
	\item{\code{var}}{The scale
	  parameter for the scale-inverse-chi-squared prior distribution
	  of the receiver effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
	  \code{receiver.var}.}
	\item{\code{var.df}}{The degrees of freedom
	  parameter for the scale-inverse-chi-squared prior distribution
	  of the receiver effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
	  \code{receiver.var.df}.}
      }
      The following parameters are associated with this term:
      \describe{
	\item{\code{receiver}}{ Numeric vector of values of each
	  vertex's random receiver effect.}
	\item{\code{receiver.var}}{ Random receiver effect's variance.}
      }
    }
    \item{\code{rsociality(var=1, var.df=3)}}{\emph{Random sociality effect.}
      Adds a random sociality effect to the model, with normal prior centered
      around \eqn{0}{0} and a variance that is estimated.
      Can be used on either a directed or an undirected network.
      The parameters are as follows:
      \describe{
	\item{\code{var}}{The scale
	  parameter for the scale-inverse-chi-squared prior distribution
	  of the sociality effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
	  \code{sociality.var}.}
	\item{\code{var.df}}{The degrees of freedom
	  parameter for the scale-inverse-chi-squared prior distribution
	  of the sociality effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
	  \code{sociality.var.df}.}
      }
      The following parameters are associated with this term:
      \describe{
	\item{\code{sociality}}{ Numeric vector of values of each
	  vertex's random sociality effect.}
	\item{\code{sociality.var}}{ Random sociality effect's variance.}
      }
    }
  }
  
  \emph{Fixed Effects}\cr
  Each coefficient for a fixed effect covariate has a normal prior whose
  mean and variance are set by the \code{mean} and \code{var} parameters
  of the term. For those formula terms that add more than one covariate,
  a vector can be given for mean and variance. If not, the vectors given
  will be repeated until the needed length is reached.

  \code{\link{ergmm}} can use model terms implemented for the
  \code{\link[=ergm-package]{ergm}} package and via the
  \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}} API. See
  \code{\link{ergm-terms}} for a list of available terms. If you wish to
  specify the prior mean and variance, you can add them to the
  call. E.g.,\cr \code{TERMNAME(..., mean=0, var=9)},\cr where
  \code{...} are the arguments for the \code{ergm} term, will initialize
  \code{TERMNAME} with prior mean of 0 and prior variance of 9.

  Some caveats:
  \itemize{
    
    \item{\code{\link[=ergm-package]{ergm}} has a binary and a valued
      mode. Regardless of the \code{\link[=families.ergmm]{family}} used,
      the \emph{binary} variant of the \code{\link[=ergm-package]{ergm}}
      term will be used in the linear predictor of the model.}
    
    \item{\code{\link[=ergm-package]{ergm}} does not support modeling
      self-loops, so terms imported in this way will always have
      predictor \code{x[i,i]==0}. This should not affect most
      situations, but if you absolutely must model self-loops and
      non-self-edges in one term, use the deprecated terms below.}
    
    \item{\code{latentnet} only fits models with dyadic
      independence. Terms that induce dyadic dependence (e.g.,
      \code{triangles}) can be used, but then the likelihood of the
      model will, effectively, be replaced with pseudolikelihood. (Note
      that under dyadic independence, the two are equal.)}
    
  }

  Each parameter in this section adds one element to \code{beta} vector.
  \describe{

    \item{\code{1(mean=0, var=9)} a.k.a. \code{intercept}
      a.k.a. \code{Intercept}}{\emph{Intercept.}  This term serves as an
      intercept term, is included by default (though, as in
      \code{\link{lm}}, it can be excluded by adding \code{+0} or
      \code{-1} into the model formula). It adds one covariate to the
      model, for which \code{x[i,j]=1} for all \code{i} and \code{j}.

      It can be used explicitly to set prior mean and variance for the
      intercept term.
	
      This term differs from the \code{ergm}'s
      \code{\link[=ergm-terms]{edges}} term if \code{g} has self-loops.}

    \item{\code{loopcov(attrname, mean=0, var=9)}}{\emph{Covariate
	effect on self-loops.}  \code{attrname} is a character string
      giving the name of a numeric (not categorical) attribute in the
      network's vertex attribute list.  This term adds one covariate
      to the model, for which \code{x[i,i]=attrname(i)} and
      \code{x[i,j]=0} for \code{i!=j}.  This term only makes sense if
      \code{g} has self-loops.}
    
    \item{\code{loopfactor(attrname, base=1, mean=0,
	var=9)}}{\emph{Factor attribute effect on self-loops.}  The
      \code{attrname} argument is a character vector giving one or
      more names of categorical attributes in the network's vertex
      attribute list. This term adds multiple covariates to the
      model, one for each of (a subset of) the unique values of the
      \code{attrname} attribute (or each combination of the
      attributes given). Each of these covariates has
      \code{x[i,i]=1} if \code{attrname(i)==l}, where \code{l} is
      that covariate's level, and \code{x[i,j]=0} otherwise. To
      include all attribute values se \code{base=0} -- because the
      sum of all such statistics equals twice the number of self-loops
      and hence a linear dependency would arise in any model also
      including \code{loops}. Thus, the \code{base} argument tells
      which value(s) (numbered in order according to the \code{sort}
      function) should be omitted. The default value, \code{base=1},
      means that the smallest (i.e., first in sorted order)
      attribute value is omitted. For example, if the \dQuote{fruit}
      factor has levels \dQuote{orange}, \dQuote{apple},
      \dQuote{banana}, and \dQuote{pear}, then to add just two
      terms, one for \dQuote{apple} and one for \dQuote{pear}, then
      set \dQuote{banana} and \dQuote{orange} to the base (remember
      to sort the values first) by using \code{nodefactor("fruit",
	base=2:3)}. For an analogous term for quantitative vertex
      attributes, see \code{nodecov}.\code{attrname} is a character
      string giving the name of a numeric (not categorical)
      attribute in the network's vertex attribute list.  This term
      adds one covariate to the model, for which
      \code{x[i,i]=attrname(i)} and \code{x[i,j]=0} for \code{i!=j}.
      This term only makes sense if \code{g} has self-loops.}

    
    \item{\code{latentcov(x, attrname=NULL, mean=0, var=9)}}{\emph{Edge covariates for the
	latent model.}

      \emph{Deprecated for networks without self-loops. Use
	\code{\link{edgecov}} instead.}
      
      \code{x} is either a matrix of
      covariates on each pair of vertices, a network, or an edge attribute on \code{g}; 
      if the latter, optional argument
      \code{attrname} provides the name of the edge attribute to
      use for edge values. \code{latentcov} can be called more
      than once, to model the effects of multiple covariates. Note that
      some covariates can be more conveniently specified using the
      following terms. 
    }
    
    \item{\code{sendercov(attrname, force.factor=FALSE, mean=0,	var=9)}}{\emph{Sender covariate effect.}

      \emph{Deprecated for networks without self-loops. Use
	\code{\link{nodeocov}}, \code{\link{nodeofactor}},
	\code{\link{nodecov}} or \code{\link{nodefactor}} instead.}

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

      \emph{Deprecated for networks without self-loops. Use
	\code{\link{nodeicov}}, \code{\link{nodeifactor}},
	\code{\link{nodecov}} or \code{\link{nodefactor}} instead.}

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

      \emph{Deprecated for networks without self-loops. Use
	\code{\link{nodecov}} instead.}
      
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
  \code{\link[ergm]{terms-ergm}}
}
\keyword{graphs}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
\keyword{nonparametric}
\keyword{cluster}
