#  File R/latentnet-package.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

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
#' produced the book cited below (DGG) that reported a comparative study of
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
#' \url{https://creativecommons.org/licenses/by-nc-nd/2.5/}.
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


#' @details
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
#' terms>}, see \code{\link{ergmTerm}}. For the list of the possible dyad
#' distribution families, see \code{\link{families.ergmm}}.
#' 
#' 
#' The arguments in the \code{\link{ergmm}} function specific to latent
#' variable models are \code{\link{ergmm.control}}. See the help page for
#' \code{\link{ergmm}} for the details.
#' 
#' The result of a latent variable model fit is an \code{\link{ergmm}} object.
#' Hence the \code{\link{summary}}, \code{print}, and \code{plot} functions
#' apply to the fits.  The \code{\link{plot.ergmm}} function has many options
#' specific to latent variable models.
#' 
#' @name latentnet-package
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
#' Statistical Software, 24(5). \doi{10.18637/jss.v024.i05}
#' 
#' Susan M. Shortreed, Mark S. Handcock, and Peter D. Hoff (2006).
#' \emph{Positional Estimation within the Latent Space Model for Networks}.
#' Methodology, 2(1), 24-33.
#' @keywords graphs package models regression nonlinear nonparametric
"_PACKAGE"





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
