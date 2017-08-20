#  File R/latentnet-package.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

#' Southern Women Data Set (Davis) as a bipartite ``network'' object
#'
#' This is a data set of 18 women observed over a nine-month period. During
#' that period, various subsets of these women had met in a series of 14
#' informal social events. The data recored which women met for which events.
#' The data is originally from Davis, Gardner and Gardner (1941) via
#' \code{UCINET} and stored as a \code{network} object.
#'
#' @name davis
#'
#' @details
#' 
#' This documentation is taken from Freeman (2003) in his usual lucid
#' description. See the reference to the paper below:
#' 
#' In the 1930s, five ethnographers, Allison Davis, Elizabeth Stubbs Davis,
#' Burleigh B. Gardner, Mary R. Gardner and J. G. St. Clair Drake, collected
#' data on stratification in Natchez, Mississippi (Warner, 1988, p. 93). They
#' produced the book cited below [DGG] that reported a comparative study of
#' social class in black and in white society. One element of this work
#' involved examining the correspondence between people's social class levels
#' and their patterns of informal interaction. DGG was concerned with the issue
#' of how much the informal contacts made by individuals were established
#' solely (or primarily) with others at approximately their own class levels.
#' To address this question the authors collected data on social events and
#' examined people's patterns of informal contacts.
#' 
#' In particular, they collected systematic data on the social activities of 18
#' women whom they observed over a nine-month period. During that period,
#' various subsets of these women had met in a series of 14 informal social
#' events. The participation of women in events was uncovered using
#' ``interviews, the records of participant observers, guest lists, and the
#' newspapers'' (DGG, p. 149). Homans (1950, p. 82), who presumably had been in
#' touch with the research team, reported that the data reflect joint
#' activities like, ``a day's work behind the counter of a store, a meeting of
#' a women's club, a church supper, a card party, a supper party, a meeting of
#' the Parent-Teacher Association, etc.''
#' 
#' This data set has several interesting properties. It is small and
#' manageable. It embodies a relatively simple structural pattern, one in
#' which, according to DGG, the women seemed to organize themselves into two
#' more or less distinct groups. Moreover, they reported that the positions -
#' core and peripheral - of the members of these groups could also be
#' determined in terms of the ways in which different women had been involved
#' in group activities. At the same time, the DGG data set is complicated
#' enough that some of the details of its patterning are less than obvious. As
#' Homans (1950, p. 84) put it, ``The pattern is frayed at the edges.'' And,
#' finally, this data set comes to us in a two-mode ``woman by event'' form.
#' Thus, it provides an opportunity to explore methods designed for direct
#' application to two-mode data. But at the same time, it can easily be
#' transformed into two one-mode matrices (woman by woman or event by event)
#' that can be examined using tools for one-mode analysis.
#' 
#' Because of these properties, this DGG data set has become something of a
#' touchstone for comparing analytic methods in social network analysis. Davis,
#' Gardner and Gardner presented an intuitive interpretation of the data, based
#' in part on their ethnographic experience in the community. Then the DGG data
#' set was picked up by Homans (1950) who provided an alternative intuitive
#' interpretation. In 1972, Phillips and Conviser used an analytic tool, based
#' on information theory, that provided a systematic way to reexamine the DGG
#' data. Since then, this data set has been analyzed again and again. It
#' reappears whenever any network analyst wants to explore the utility of some
#' new tool for analyzing data.
#' 
#' If the source of the data set does not specified otherwise, this data set is
#' protected by the Creative Commons License
#' \url{http://creativecommons.org/licenses/by-nc-nd/2.5/}.
#' 
#' When publishing results obtained using this data set the original authors
#' should be cited.  In addition this package should be cited.
#' 
#' @seealso statnet, network, ergm, ergm
#' @references Davis, A., Gardner, B. B. and M. R. Gardner (1941) \emph{Deep
#' South,} Chicago: The University of Chicago Press.
#' 
#' Linton C. Freeman (2003). \emph{Finding Social Groups: A Meta-Analysis of
#' the Southern Women Data}, In Ronald Breiger, Kathleen Carley and Philippa
#' Pattison, eds. Dynamic Social Network Modeling and Analysis. Washington: The
#' National Academies Press.
#' @source Linton C. Freeman (2003).  \emph{Finding Social Groups: A
#' Meta-Analysis of the Southern Women Data}, In Ronald Breiger, Kathleen
#' Carley and Philippa Pattison, eds. Dynamic Social Network Modeling and
#' Analysis. Washington: The National Academies Press.
#' @keywords data
#' @docType data
#' @examples
#' 
#' \donttest{
#' data(davis)
#' # Fit a 2D 2-cluster fit and plot.
#' davis.fit<-ergmm(davis~euclidean(d=2,G=2)+rsociality)
#' plot(davis.fit,pie=TRUE,rand.eff="sociality")
#' }
#' 
NULL




#' Class of Fitted Exponential Random Graph Mixed Models
#'
#' A class \code{\link[=ergmm.object]{ergmm}} to represent a fitted exponential
#' random graph mixed model. The output of \code{\link{ergmm}}.
#' 
#'  There are methods \code{\link{summary.ergmm}}, \code{print.ergmm},
#'  \code{\link{plot.ergmm}}, \code{\link{predict.ergmm}}, and
#'  \code{\link{as.mcmc.list.ergmm}}.
#'  
#'  The structure of \code{ergmm} is as follows:
#'  \describe{
#'     \item{\code{sample}}{ An object of class \code{\link[=ergmm.par.list.object]{ergmm.par.list}} containing the
#'    MCMC sample from the posterior. If the run had multiple threads, their output is concatenated.}
#'     \item{\code{mcmc.mle}}{ A list containing the parameter
#'    configuration of the highest-likelihood MCMC iteration. }
#'     \item{\code{mcmc.pmode}}{ A list containing the parameter
#'    configuration of the highest-joint-density (conditional on cluster
#'    assignments) MCMC iteration. }
#'     \item{\code{mkl}}{ A list containing the MKL estimate. }
#'     \item{\code{model}}{ A list containing the model
#'    that was fitted.}
#'     \item{\code{prior}}{ A list containing the
#'    information about the prior distribution used. It can be passed as
#'    parameter \code{prior} to \code{\link{ergmm}} to reproduce the prior
#'    in a new fit.}
#'     \item{\code{control}}{ A list containing the
#'    information about the model fit settings that do not affect the
#'    posterior distribution. It can be passed as
#'    parameter \code{control} to \code{\link{ergmm}} to reproduce control
#'    parameters in a new fit.}
#'     \item{\code{mle}}{ A list containing the MLE,
#'    conditioned on cluster assignments.}
#'     \item{\code{pmode}}{ A list containing the posterior mode,
#'    conditioned on cluster assignments.}
#'     \item{\code{burnin.start}}{ A list containing the starting
#'    value for the burnin.}
#'     \item{\code{main.start}}{  A list (or a list of lists, for a
#'    multithreaded run) containing the starting
#'    value for the sampling. }
#'  }
#' 
#' @name ergmm-class
#' @aliases ergmm.object print.ergmm show.ergmm
#' @seealso \code{\link{ergmm}}, \code{\link{summary.ergmm}},
#' \code{\link{plot.ergmm}}, \code{\link{predict.ergmm}},
#' \code{\link{as.mcmc.list.ergmm}}
#' @keywords graphs regression models
NULL








#' Edge Weight Distribution Families
#'
#' @name ergmm-families
#'
#' @description
#' Family-link combinations supported by \code{\link{ergmm}}.
#'
#' @details
#' Each supported family has a family of functions, of the form \code{pY.}-,
#' \code{lpY.}-, \code{EY.}-, \code{dlpY.deta.}-, \code{dlpY.ddispersion.}-,
#' \code{lpYc.}-, \code{rsm.}-, followed by the family's name, for the
#' respective family's name, representing the family's likelihood,
#' log-likelihood, expectation, derivative of log-likelihood with repect to the
#' linear predictor, derivative of log-likelihood with respect to the
#' dispersion parameter, log-normalizing-constant, and random sociomatrix
#' generation functions.
#' 
#' On the \code{C} side, similar functions exist, but becuase of static typing,
#' are also provided for ``continuous'' versions of those families. These
#' should not be used on their own, but are used in estimating MKL positions
#' from the posterior distribution.
#' 
#' @aliases families.ergmm family ergmm.families dlpY.deta.Bernoulli.logit
#' dlpY.deta.binomial.logit dlpY.deta.fs dlpY.deta.Poisson.log
#' dlpY.deta.normal.identity dlpY.ddispersion.fs
#' dlpY.ddispersion.normal.identity lpYc.Bernoulli.logit lpYc.binomial.logit
#' lpYc.normal.identity lpYc.fs lpYc.Poisson.log lpY.Bernoulli.logit
#' lpY.binomial.logit lpY.Poisson.log lpY.normal.identity lpY.fs
#' EY.Bernoulli.logit EY.binomial.logit EY.fs EY.Poisson.log EY.normal.identity
#' pY.fs pY.Poisson.log pY.Bernoulli.logit pY.binomial.logit pY.normal.identity
#' rsm.fs rsm.Poisson.log rsm.binomial.logit rsm.Bernoulli.logit
#' rsm.normal.identity family.IDs family.names fam.par.check
#' 
#' @section Family-link combinations:
#'
#'   \tabular{rlllll}{
#'    ID \tab \code{C} name          \tab \R name   \tab Type       \tab Family    \tab Link  \cr
#'     1 \tab \code{Bernoulli_logit}       \tab \code{Bernoulli.logit} \tab Discrete   \tab Bernoulli \tab logit \cr
#'     2 \tab \code{binomial_logit}        \tab \code{binomial.logit}  \tab Discrete   \tab binomial  \tab logit \cr
#'     3 \tab \code{Poisson_log}           \tab \code{Poisson.log}     \tab Discrete   \tab Possion   \tab log   \cr
#'     4 \tab \code{Bernoulli_cont_logit} \tab NA              \tab Continuous \tab Bernoulli \tab logit \cr
#'     5 \tab \code{binomial_cont_logit}  \tab NA              \tab Continuous \tab binomial  \tab logit \cr
#'     6 \tab \code{Poisson_cont_log}     \tab NA              \tab Continuous \tab Possion   \tab log \cr  
#'     7 \tab \code{normal_identity}      \tab \code{normal.identity} \tab Continuous \tab normal   \tab identity  
#'   }
#'   \code{.link} can be omited when not ambiguous. Some families
#'   require an appropriate \code{fam.par} argument to be supplied to
#'   \code{\link{ergmm}}:
#'   \describe{
#'     \item{binomial families}{a mandatory \code{trials} parameter for the
#'       number of trials (same for every dyad) whose success the response
#'       counts represent}
#'     \item{normal}{a mandatory \code{prior.var} and \code{prior.var.df} parameter for the prior scale and degrees of freedom of the variance of
#'       the dyad values}
#'   }
#' @keywords graphs models regression
NULL


## #' Internal latentnet Objects
## #' 
## #' Internal latentnet functions.
## #' 
## #' 
## #' @aliases mbc.VII.EM print.ergmm.model ergmm.model ergmm.model.object getYm
## #' bipartite.augment observed.dyads statsreeval.ergmm ergmm.lpY ergmm.lpY.C
## #' ergmm.lpY.grad ergmm.lp ergmm.lpBeta ergmm.lpBeta.grad ergmm.lp.grad
## #' ergmm.lp.grad.approx ergmm.lpLV ergmm.lpLV.grad ergmm.lpRE ergmm.lpRE.grad
## #' ergmm.lpREV ergmm.lpREV.grad ergmm.lpZ ergmm.lpZ.grad ergmm.lpdispersion
## #' ergmm.lpdispersion.grad ergmm.eta ergmm.EY lp.works lpsum cov.beta.ext
## #' get.init.deltas get.sample.deltas find.mle find.mle.loop find.mpe
## #' find.pmode.loop find.mkl add.mcmc.mle.mle.ergmm add.mcmc.pmode.pmode.ergmm
## #' add.mkl.mbc.ergmm add.mkl.pos.ergmm xtabs.ergmm clust.homogeneity
## #' @seealso ergmm
## #' @keywords internal
## NULL




#' Latent position and cluster models for networks
#'
#'
#' The package \code{latentnet} is used to fit latent cluster random effect
#' models, where the probability of a network \eqn{g}, on a set of nodes is a
#' product of dyad probabilities, each of which is a GLM with linear component
#' \eqn{\eta_{i,j}=\sum_{k=1}^p \beta_k
#' X_{i,j,k}+d(Z_i,Z_j)+\delta_i+\gamma_j}, where \eqn{X} is an array of dyad
#' covariates, \eqn{\beta} is a vector of covariate coefficients, \eqn{Z_i} is
#' the latent space position of node \eqn{i}, \eqn{d(\cdot,\cdot)} is a
#' function of the two positions: either negative Euclidean
#' (\eqn{-||Z_i-Z_j||}) or bilinear (\eqn{Z_i\cdot Z_j}), and \eqn{\delta} and
#' \eqn{\gamma} are vectors of sender and receiver effects. (Note that these
#' are different from the eigenmodel of Hoff (2007) ``Modeling homophily and
#' stochastic equivalence in symmetric relational data'', fit by package
#' \code{eigenmodel}.)
#' 
#' The \code{\link{ergmm}} specifies models via: \code{g ~ <model terms>} where
#' \code{g} is a \code{network} object For the list of possible \code{<model
#' terms>}, see \link{terms.ergmm}. For the list of the possible dyad
#' distribution families, see \link{families.ergmm}.
#' 
#' 
#' The arguments in the \code{\link{ergmm}} function specific to latent
#' variable models are \code{ergmm.control}. See the help page for
#' \code{\link{ergmm}} for the details.
#' 
#' The result of a latent variable model fit is an \code{\link{ergmm}} object.
#' Hence the \code{\link{summary}}, \code{print}, and \code{plot} functions
#' apply to the fits.  The \code{\link{plot.ergmm}} function has many options
#' specific to latent variable models. % are \code{mle}, \code{pie},
#' \code{contour}, % \code{contour.colors}.  See the help page for
#' \code{\link{plot.ergmm}} for the details.
#' 
#' @name latentnet-package
#' @aliases latentnet-package latentnet
#' @return \code{\link{ergmm}} returns an object of class 'ergmm' that is a
#' list.
#' @seealso \code{\link{ergmm}}, \link{terms.ergmm}
#' @references Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum (2007).
#' \emph{Model-Based Clustering for Social Networks}.  Journal of the Royal
#' Statistical Society: Series A (Statistics in Society), 170(2), 301-354.
#' 
#' Peter D. Hoff (2005). \emph{Bilinear Mixed Effects Models for Dyadic Data}.
#' Journal of the American Statistical Association, 100(469), 286-295.
#' 
#' Peter D. Hoff, Adrian E. Raftery and Mark S. Handcock (2002).  \emph{Latent
#' space approaches to social network analysis}.  Journal of the American
#' Statistical Association, 97(460), 1090-1098.
#' 
#' Pavel N. Krivitsky, Mark S. Handcock, Adrian E. Raftery, and Peter D. Hoff
#' (2009).  \emph{Representing degree distributions, clustering, and homophily
#' in social networks with latent cluster random effects models}.  Social
#' Networks, 31(3), 204-213.
#' 
#' Pavel N. Krivitsky and Mark S. Handcock (2008).  \emph{Fitting Position
#' Latent Cluster Models for Social Networks with \code{latentnet}}. Journal of
#' Statistical Software, 24(5).
#' 
#' Susan M. Shortreed, Mark S. Handcock, and Peter D. Hoff (2006).
#' \emph{Positional Estimation within the Latent Space Model for Networks}.
#' Methodology, 2(1), 24-33.
#' @docType package
#' @keywords graphs package models regression nonlinear nonparametric
NULL







#' Model Terms for Latent Space Random Graph Model
#'
#' Model terms that can be used in an \code{\link{ergmm}} formula and their
#' parameter names.
#' 
#' @name ergmm-terms
#'
#' @section Model Terms:
#' 
#'   The \code{\link{latentnet}} package itself allows only
#'  dyad-independent terms. In the formula for the model, the model terms are various function-like
#'  calls, some of which require arguments, separated by \code{+} signs.
#'
#'  \emph{Latent Space Effects}\cr
#'  \describe{
#'    \item{\code{euclidean(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
#'  mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)}}{\emph{(Negative)
#'  Euclidean distance model term, with
#'	optional clustering.}
#'      Adds a term to the model equal to the negative Eucledean distance
#'      \eqn{-||Z_i-Z_j||}{-dist(Z[i],Z[j])}, where \eqn{Z_i}{Z[i]} and \eqn{Z_j}{Z[j]}
#'      are the positions of their respective actors in an unobserved social
#'      space. These positions may optionally have a finite spherical
#'      Gaussian mixture clustering structure. This term was previously
#'      called \code{latent} which now fits negative Euclidean latent
#'      space model with a warning.
#'      The parameters are as follows:
#'      \describe{
#'	\item{\code{d}}{The dimension of the latent space.}
#'	\item{\code{G}}{The number of groups (0 for no clustering).}
#'	\item{\code{var.mul}}{In the absence of \code{var}, this
#'	  argument will be used as a scaling factor for a
#'	  function of average cluster size and latent space dimension to
#'	  set \code{var}. To set it in the \code{prior} argument to \code{\link{ergmm}}, use \code{Z.var.mul}.}
#'	\item{\code{var}}{If given, the scale
#'	  parameter for the scale-inverse-chi-squared prior distribution of the within-cluster
#'	  variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{Z.var}.}
#'	\item{\code{var.df.mul}}{In the absence of \code{var.df}, this
#'	  argument is the multiplier for the square root of average
#'	  cluster size, which serves in place of \code{var.df}. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{Z.var.df.mul}.}
#'	\item{\code{var.df}}{The degrees of freedom
#'	  parameter for the scale-inverse-chi-squared prior distribution of the within-cluster
#'	  variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{Z.var.df}.}
#'	\item{\code{mean.var.mul}}{In the absence of \code{mean.var},
#'	  the multiplier for a function of number of vertices and latent space
#'	  dimension to set \code{mean.var}. To set it in the
#'	  \code{prior} argument to \code{\link{ergmm}}, use \code{Z.mean.var.mul}.}
#'	\item{\code{mean.var}}{The variance of
#'	  the spherical Gaussian prior distribution of the cluster means. To set it in the
#'	  \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{Z.mean.var}.}
#'	\item{\code{pK.mul}}{In the absence of \code{pK}, this argument
#'	  is the multiplier for the square root of the average cluster size,
#'	  which is used as \code{pK}. To set it in
#'	  the \code{prior} argument to \code{\link{ergmm}}, use \code{Z.pK}.}
#'	\item{\code{pK}}{The parameter of the Dirichilet prior
#'	  distribution of cluster assignment probabilities. To set it in
#'	  the \code{prior} argument to \code{\link{ergmm}}, use \code{Z.pK}.}
#'      }
#'      The following parameters are associated with this term:
#'      \describe{
#'	\item{\code{Z}}{ Numeric matrix with rows being latent space
#'	  positions.}
#'	\item{\code{Z.K} (when \eqn{\code{G}>0})}{ Integer vector of cluster assignments. }
#'	\item{\code{Z.mean} (when \eqn{\code{G}>0})}{ Numeric matrix with rows being cluster means. }
#'	\item{\code{Z.var}  (when \eqn{\code{G}>0})}{ Depending on the model, either a numeric vector with
#'	  within-cluster variances or a numeric scalar with the overal latent space variance. }
#'	\item{\code{Z.pK}  (when \eqn{\code{G}>0})}{ Numeric vector of probabilities of a vertex being
#'	  in a particular cluster.} 
#'      }
#'    }
#'    \item{\code{bilinear(d, G=0, var.mul=1/8, var=NULL, var.df.mul=1, var.df=NULL,
#'  mean.var.mul=1, mean.var=NULL, pK.mul=1, pK=NULL)}}{\emph{Bilinear
#'  latent model term, with
#'	optional clustering.}
#'      Adds a term to the model equal to the inner product of the latent positions:
#'      \eqn{Z_i \cdot Z_j}{sum(Z[i]*Z[j])}, where \eqn{Z_i}{Z[i]} and \eqn{Z_j}{Z[j]}
#'      are the positions of their respective actors in an unobserved social
#'      space. These positions may optionally have a finite spherical
#'      Gaussian mixture clustering structure. \emph{Note: For a
#'	bilinear latent space effect, two actors being closer in the
#'	clustering sense does not necessarily mean that the expected value of
#'	a tie between them is higher. Thus, a warning is printed when
#'	this model is combined with clustering.}
#'      The parameters are as follows:
#'      \describe{
#'	\item{\code{d}}{The dimension of the latent space.}
#'	\item{\code{G}}{The number of groups (0 for no clustering).}
#'	\item{\code{var.mul}}{In the absence of \code{var}, this
#'	  argument will be used as a scaling factor for a
#'	  function of average cluster size and latent space dimension to
#'	  set \code{var}. To set it in the \code{prior} argument to \code{\link{ergmm}}, use \code{Z.var.mul}.}
#'	\item{\code{var}}{If given, the scale
#'	  parameter for the scale-inverse-chi-squared prior distribution of the within-cluster
#'	  variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{Z.var}.}
#'	\item{\code{var.df.mul}}{In the absence of \code{var.df}, this
#'	  argument is the multiplier for the square root of average
#'	  cluster size, which serves in place of \code{var.df}. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{Z.var.df.mul}.}
#'	\item{\code{var.df}}{The degrees of freedom
#'	  parameter for the scale-inverse-chi-squared prior distribution of the within-cluster
#'	  variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{Z.var.df}.}
#'	\item{\code{mean.var.mul}}{In the absence of \code{mean.var},
#'	  the multiplier for a function of number of vertices and latent space
#'	  dimension to set \code{mean.var}. To set it in the
#'	  \code{prior} argument to \code{\link{ergmm}}, use \code{Z.mean.var.mul}.}
#'	\item{\code{mean.var}}{The variance of
#'	  the spherical Gaussian prior distribution of the cluster means. To set it in the
#'	  \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{Z.mean.var}.}
#'	\item{\code{pK.mul}}{In the absence of \code{pK}, this argument
#'	  is the multiplier for the square root of the average cluster size,
#'	  which is used as \code{pK}. To set it in
#'	  the \code{prior} argument to \code{\link{ergmm}}, use \code{Z.pK}.}
#'	\item{\code{pK}}{The parameter of the Dirichilet prior
#'	  distribution of cluster assignment probabilities. To set it in
#'	  the \code{prior} argument to \code{\link{ergmm}}, use \code{Z.pK}.}
#'      }
#'      The following parameters are associated with this term:
#'      \describe{
#'	\item{\code{Z}}{ Numeric matrix with rows being latent space
#'	  positions.}
#'	\item{\code{Z.K} (when \eqn{\code{G}>0})}{ Integer vector of cluster assignments. }
#'	\item{\code{Z.mean} (when \eqn{\code{G}>0})}{ Numeric matrix with rows being cluster means. }
#'	\item{\code{Z.var}  (when \eqn{\code{G}>0})}{ Depending on the model, either a numeric vector with
#'	  within-cluster variances or a numeric scalar with the overal latent space variance. }
#'	\item{\code{Z.pK}  (when \eqn{\code{G}>0})}{ Numeric vector of probabilities of a vertex being
#'	  in a particular cluster.} 
#'      }
#'    }
#'  }
#'  \emph{Actor-specific effects}\cr
#'  \describe{
#'    \item{\code{rsender(var=1, var.df=3)}}{\emph{Random sender effect.}
#'      Adds a random sender effect to the model, with normal prior centered
#'      around \eqn{0}{0} and a variance that is estimated.
#'      Can only be used on directed networks.
#'      The parameters are as follows:
#'      \describe{
#'	\item{\code{var}}{The scale
#'	  parameter for the scale-inverse-chi-squared prior distribution
#'	  of the sender effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{sender.var}.}
#'	\item{\code{var.df}}{The degrees of freedom
#'	  parameter for the scale-inverse-chi-squared prior distribution
#'	  of the sender effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{sender.var.df}.
#'	}
#'      }
#'      The following parameters are associated with this term:
#'      \describe{
#'	\item{\code{sender}}{ Numeric vector of values of each
#'	  vertex's random sender effect.}
#'	\item{\code{sender.var}}{ Random sender effect's variance.}
#'      }
#'    }
#'    \item{\code{rreceiver(var=1, var.df=3)}}{\emph{Random receiver effect.}
#'      Adds a random receiver effect to the model, with normal prior centered
#'      around \eqn{0}{0} and a variance that is estimated.
#'      Can only be used on directed networks.
#'      The parameters are as follows:
#'      \describe{
#'	\item{\code{var}}{The scale
#'	  parameter for the scale-inverse-chi-squared prior distribution
#'	  of the receiver effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{receiver.var}.}
#'	\item{\code{var.df}}{The degrees of freedom
#'	  parameter for the scale-inverse-chi-squared prior distribution
#'	  of the receiver effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{receiver.var.df}.}
#'      }
#'      The following parameters are associated with this term:
#'      \describe{
#'	\item{\code{receiver}}{ Numeric vector of values of each
#'	  vertex's random receiver effect.}
#'	\item{\code{receiver.var}}{ Random receiver effect's variance.}
#'      }
#'    }
#'    \item{\code{rsociality(var=1, var.df=3)}}{\emph{Random sociality effect.}
#'      Adds a random sociality effect to the model, with normal prior centered
#'      around \eqn{0}{0} and a variance that is estimated.
#'      Can be used on either a directed or an undirected network.
#'      The parameters are as follows:
#'      \describe{
#'	\item{\code{var}}{The scale
#'	  parameter for the scale-inverse-chi-squared prior distribution
#'	  of the sociality effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{sociality.var}.}
#'	\item{\code{var.df}}{The degrees of freedom
#'	  parameter for the scale-inverse-chi-squared prior distribution
#'	  of the sociality effect variance. To set it in the \code{prior} argument to \code{\link{ergmm}}, use
#'	  \code{sociality.var.df}.}
#'      }
#'      The following parameters are associated with this term:
#'      \describe{
#'	\item{\code{sociality}}{ Numeric vector of values of each
#'	  vertex's random sociality effect.}
#'	\item{\code{sociality.var}}{ Random sociality effect's variance.}
#'      }
#'    }
#'  }
#'  
#'  \emph{Fixed Effects}\cr
#'  Each coefficient for a fixed effect covariate has a normal prior whose
#'  mean and variance are set by the \code{mean} and \code{var} parameters
#'  of the term. For those formula terms that add more than one covariate,
#'  a vector can be given for mean and variance. If not, the vectors given
#'  will be repeated until the needed length is reached.
#'
#'  \code{\link{ergmm}} can use model terms implemented for the
#'  \code{\link[=ergm-package]{ergm}} package and via the
#'  \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}} API. See
#'  \code{\link{ergm-terms}} for a list of available terms. If you wish to
#'  specify the prior mean and variance, you can add them to the
#'  call. E.g.,\cr \code{TERMNAME(..., mean=0, var=9)},\cr where
#'  \code{...} are the arguments for the \code{ergm} term, will initialize
#'  \code{TERMNAME} with prior mean of 0 and prior variance of 9.
#'
#'  Some caveats:
#'  \itemize{
#'    
#'    \item{\code{\link[=ergm-package]{ergm}} has a binary and a valued
#'      mode. Regardless of the \code{\link[=families.ergmm]{family}} used,
#'      the \emph{binary} variant of the \code{\link[=ergm-package]{ergm}}
#'      term will be used in the linear predictor of the model.}
#'    
#'    \item{\code{\link[=ergm-package]{ergm}} does not support modeling
#'      self-loops, so terms imported in this way will always have
#'      predictor \code{x[i,i]==0}. This should not affect most
#'      situations, but if you absolutely must model self-loops and
#'      non-self-edges in one term, use the deprecated terms below.}
#'    
#'    \item{\code{latentnet} only fits models with dyadic
#'      independence. Terms that induce dyadic dependence (e.g.,
#'      \code{triangles}) can be used, but then the likelihood of the
#'      model will, effectively, be replaced with pseudolikelihood. (Note
#'      that under dyadic independence, the two are equal.)}
#'    
#'  }
#'
#'  Each parameter in this section adds one element to \code{beta} vector.
#'  \describe{
#'
#'    \item{\code{1(mean=0, var=9)} a.k.a. \code{intercept}
#'      a.k.a. \code{Intercept}}{\emph{Intercept.}  This term serves as an
#'      intercept term, is included by default (though, as in
#'      \code{\link{lm}}, it can be excluded by adding \code{+0} or
#'      \code{-1} into the model formula). It adds one covariate to the
#'      model, for which \code{x[i,j]=1} for all \code{i} and \code{j}.
#'
#'      It can be used explicitly to set prior mean and variance for the
#'      intercept term.
#'	
#'      This term differs from the \code{ergm}'s
#'      \code{\link[=ergm-terms]{edges}} term if \code{g} has self-loops.}
#'
#'    \item{\code{loopcov(attrname, mean=0, var=9)}}{\emph{Covariate
#'	effect on self-loops.}  \code{attrname} is a character string
#'      giving the name of a numeric (not categorical) attribute in the
#'      network's vertex attribute list.  This term adds one covariate
#'      to the model, for which \code{x[i,i]=attrname(i)} and
#'      \code{x[i,j]=0} for \code{i!=j}.  This term only makes sense if
#'      \code{g} has self-loops.}
#'    
#'    \item{\code{loopfactor(attrname, base=1, mean=0,
#'	var=9)}}{\emph{Factor attribute effect on self-loops.}  The
#'      \code{attrname} argument is a character vector giving one or
#'      more names of categorical attributes in the network's vertex
#'      attribute list. This term adds multiple covariates to the
#'      model, one for each of (a subset of) the unique values of the
#'      \code{attrname} attribute (or each combination of the
#'      attributes given). Each of these covariates has
#'      \code{x[i,i]=1} if \code{attrname(i)==l}, where \code{l} is
#'      that covariate's level, and \code{x[i,j]=0} otherwise. To
#'      include all attribute values se \code{base=0} -- because the
#'      sum of all such statistics equals twice the number of self-loops
#'      and hence a linear dependency would arise in any model also
#'      including \code{loops}. Thus, the \code{base} argument tells
#'      which value(s) (numbered in order according to the \code{sort}
#'      function) should be omitted. The default value, \code{base=1},
#'      means that the smallest (i.e., first in sorted order)
#'      attribute value is omitted. For example, if the \dQuote{fruit}
#'      factor has levels \dQuote{orange}, \dQuote{apple},
#'      \dQuote{banana}, and \dQuote{pear}, then to add just two
#'      terms, one for \dQuote{apple} and one for \dQuote{pear}, then
#'      set \dQuote{banana} and \dQuote{orange} to the base (remember
#'      to sort the values first) by using \code{nodefactor("fruit",
#'	base=2:3)}. For an analogous term for quantitative vertex
#'      attributes, see \code{nodecov}.\code{attrname} is a character
#'      string giving the name of a numeric (not categorical)
#'      attribute in the network's vertex attribute list.  This term
#'      adds one covariate to the model, for which
#'      \code{x[i,i]=attrname(i)} and \code{x[i,j]=0} for \code{i!=j}.
#'      This term only makes sense if \code{g} has self-loops.}
#'
#'    
#'    \item{\code{latentcov(x, attrname=NULL, mean=0, var=9)}}{\emph{Edge covariates for the
#'	latent model.}
#'
#'      \emph{Deprecated for networks without self-loops. Use
#'	\code{\link{edgecov}} instead.}
#'      
#'      \code{x} is either a matrix of
#'      covariates on each pair of vertices, a network, or an edge attribute on \code{g}; 
#'      if the latter, optional argument
#'      \code{attrname} provides the name of the edge attribute to
#'      use for edge values. \code{latentcov} can be called more
#'      than once, to model the effects of multiple covariates. Note that
#'      some covariates can be more conveniently specified using the
#'      following terms. 
#'    }
#'    
#'    \item{\code{sendercov(attrname, force.factor=FALSE, mean=0,	var=9)}}{\emph{Sender covariate effect.}
#'
#'      \emph{Deprecated for networks without self-loops. Use
#'	\code{\link{nodeocov}}, \code{\link{nodeofactor}},
#'	\code{\link{nodecov}} or \code{\link{nodefactor}} instead.}
#'
#'      \code{attrname} is a character string giving the name of an
#'      attribute in the network's vertex attribute list.
#'      If the attribute is numeric, This term adds one covariate
#'      to the model equaling \code{attrname(i)}. If the attribute is not
#'      numeric or \code{force.factor==TRUE}, this term adds \eqn{p-1}
#'      covariates to the model,
#'      where \eqn{p} is the number of unique values of \code{attrname}.
#'      The \eqn{k}th such covariate has the value \code{attrname(i) == value(k+1)}, where
#'      \code{value(k)} is the \eqn{k}th smallest unique value of the
#'      \code{attrname} attribute. This term only makes
#'      sense if \code{g} is directed.}
#'    
#'    \item{\code{receivercov(attrname, force.factor=FALSE, mean=0, var=9)}}{\emph{Receiver covariate effect.}
#'
#'      \emph{Deprecated for networks without self-loops. Use
#'	\code{\link{nodeicov}}, \code{\link{nodeifactor}},
#'	\code{\link{nodecov}} or \code{\link{nodefactor}} instead.}
#'
#'      \code{attrname} is a character string giving the name of an
#'      attribute in the network's vertex attribute list.
#'      If the attribute is numeric, This term adds one covariate
#'      to the model equaling \code{attrname(j)}. If the attribute is not
#'      numeric or \code{force.factor==TRUE}, this term adds \eqn{p-1}
#'      covariates to the model,
#'      where \eqn{p} is the number of unique values of \code{attrname}.
#'      The \eqn{k}th such covariate has the value \code{attrname(j) == value(k+1)}, where
#'      \code{value(k)} is the \eqn{k}th smallest unique value of the
#'      \code{attrname} attribute. This term only makes
#'      sense if \code{g} is directed.}
#'    
#'    \item{\code{socialitycov(attrname, force.factor=FALSE, mean=0, var=9)}}{\emph{Sociality covariate effect.}
#'
#'      \emph{Deprecated for networks without self-loops. Use
#'	\code{\link{nodecov}} instead.}
#'      
#'      \code{attrname} is a character string giving the name of an
#'      attribute in the network's vertex attribute list.
#'      If the attribute is numeric, This term adds one covariate
#'      to the model equaling \code{attrname(i)+attrname(j)}. If the attribute is not
#'      numeric or \code{force.factor==TRUE}, this term adds \eqn{p-1}
#'      covariates to the model,
#'      where \eqn{p} is the number of unique values of \code{attrname}.
#'      The \eqn{k}th such covariate has the value \code{attrname(i) ==
#'	value(k+1) + attrname(j) == value(k+1)}, where
#'      \code{value(k)} is the \eqn{k}th smallest unique value of the
#'      \code{attrname} attribute. This term makes sense whether or not
#'      \code{g} is directed.}
#'  }
#'
#' 
#' @aliases terms.ergmm ergmm.terms terms-ergmm ergmm-terms dlpY.dZ.bilinear
#' dlpY.dZ.negative.Euclidean dlpY.dZ.negative.Euclidean2 dlpY.dZ.fs
#' latent.effect.fs latent.effect.IDs latent.effect.invariances
#' latent.effect.names InitErgmm.intercept InitErgmm.Intercept InitErgmm.1 InitErgmm.latent
#' InitErgmm.loops InitErgmm.loopcov InitErgmm.loopfactor InitErgmm.euclidean
#' InitErgmm.bilinear InitErgmm.euclidean2 InitErgmm.latentcov
#' InitErgmm.receivercov InitErgmm.rreceiver InitErgmm.rsender
#' InitErgmm.rsociality InitErgmm.sendercov InitErgmm.socialitycov intercept
#' Intercept 1 loops loopcov loopfactor euclidean bilinear euclidean2 latentcov
#' receivercov rreceiver rsender rsociality sendercov socialitycov
#' latent.effect.bilinear latent.effect.negative.Euclidean
#' latent.effect.negative.Euclidean2
#'
#' @section Model Terms:
#'
#'
#'
#'
#' @seealso \code{\link{ergmm}} \code{\link[ergm]{terms-ergm}}
#' @keywords graphs models regression nonlinear nonparametric cluster
NULL





#' Read Highland Tribes
#' 
#' A network of political alliances and enmities among the 16 Gahuku-Gama
#' sub-tribes of Eastern Central Highlands of New Guinea, documented by Read
#' (1954).
#' 
#' This network shows 3 clusters.
#'
#' @name tribes
#' @docType data
#' @format
#'  An undirected \code{\link[network]{network}} object with no loops, having the following attributes:
#'  \describe{
#'    \item{\code{\%v\% "vertex.names"}}{Character attribute with names of tribes.}
#'    \item{\code{\%e\% "pos"}}{Logical attribute indicating an
#'      alliance relationship.}
#'    \item{\code{\%e\% "neg"}}{Logical attribute indicating a hostile
#'      relationship ("rova").}
#'    \item{\code{\%e\% "sign"}}{Numeric attribute coding -1 for enmity, 0
#'      for no relationship, and 1 for alliance.}
#'    \item{\code{\%e\% "sign.012"}}{Numeric attribute coding 0 for enmity, 1
#'      for no relationship, and 2 for alliance.}
#'  }
#'  Because of limitations of \code{\link[network]{network}} objects, the object
#'  itself is a complete graph, and is thus meaningless if used directly
#'  or plotted.
#' @references Taken from UCINET IV, which cites the following: Hage P. and
#' Harary F. (1983). Structural models in anthropology. Cambridge: Cambridge
#' University Press. (See p 56-60).  Read K. (1954). Cultures of the central
#' highlands, New Guinea. Southwestern Journal of Anthropology, 10, 1-43.
#' @source
#' \url{http://vlado.fmf.uni-lj.si/pub/networks/data/UciNet/UciData.htm#gama},
#' with corrections from Read (1954).
#' @keywords multivariate cluster graphs
#' @examples
#' 
#' \donttest{
#' data(tribes)
#' # Only model positive ties:
#' tribes.fit<-ergmm(tribes~euclidean(d=2,G=3),response="pos")
#' # Edge color must be set manually, for green ties to represent alliance
#' # and for red ties to represent enmity.
#' plot(tribes.fit,edge.col=as.matrix(tribes,"pos",m="a")*3+as.matrix(tribes,"neg",m="a")*2,pie=TRUE)
#' # Model both positive and negative ties:
#' tribes.fit3<-ergmm(tribes~euclidean(d=2,G=3),response="sign.012",
#'                    family="binomial.logit",fam.par=list(trials=2))
#' # Edge color must be set manually, for green ties to represent alliance
#' # and for red ties to represent enmity.
#' plot(tribes.fit3,edge.col=as.matrix(tribes,"pos",m="a")*3+as.matrix(tribes,"neg",m="a")*2,pie=TRUE)
#' }
#' 
NULL
