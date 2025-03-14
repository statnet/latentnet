# latentnet 2.12.0

* Binomial family can now take a sociomatrix for the number of trials, allowing a different number for each dyad.

# latentnet 2.11.0

* Term lookup is improved, using `statnet.common`'s function locator, cleaning up the package name space.
* Integration of `latentnet`'s term documentation into `ergm`'s.
* Updates to work better with the latest version of `Roxygen`.
* Documentation issues flagged by CRAN fixed.

# latentnet 2.10.6

* Miscellaneous fixes and updates for compatibility to ergm 4.2.
* Miscellaneous documentation improvements.

# latentnet 2.10.5

* Fix for a C-level bug exposed by link-time optimization.
* Fix for a C-level memory issue exposed by UBSAN.
* Miscelanious fixes for code analysis warnings.
* Miscellaneous documentation improvements.

# latentnet 2.10.1

* Changes for compatibility with GCC 10.
* latentnet can now fit models with no fixed effects.

# latentnet 2.9.0

* Changes for compatibility with statnet.common 3.1 and ergm 3.9.4.

# latentnet 2.8.0

* A bug in the MCMC sampling code, exposed by Philip Leifield's work and localized by Mark Handcock and Peter Hoff, has been fixed. The bug appears to have been introduced in latentnet 2.0 or earlier.
* This bug affects models for directed networks which have a fixed effect that is asymmetrical. In particular, uses of nodeocov(), nodeicov(), edgecov(), and latentcov() with an asymmetrical covariate matrix are affected. The bug has the effect of transposing the covariate matrix: nodeocov() effectively becomes nodeicov() and vice versa, and edgecov()'s matrix is transposed. See help("ergmm-terms") for more information.
* Term latent() has been added back, producing a sensible error message.
* The package now uses Roxygen to manage namespaces and documentation.
* Some miscellaneous bugs were also fixed; in particular, by Jake Fisher and Jordan T. Bates.

# latentnet 2.7.2

* A bug in model initialization, reported by Ryan Haunfelder, has been fixed. It appears to have been introduced in version 2.5.1.  Only undirected networks are affected, if the model uses any of the following terms:
    * latentcov()
    * socialitycov()
    * loopfactor()
    * loopcov()
    * intercept specified implicitly or via 1 in the formula (though not one specified via edges term)
* The covariate (X) values associated with the affected term was doubled. This means that (for a weak prior), the coefficients and the standard errors for the affected terms were both halved, with other terms, including latent ones, unaffected. Significance and fit metrics were not affected, and neither were predictions.
* Behavior of latentcov() for an undirected network given an asymmetric covariate matrix has changed: where before, the matrix were symmetrized by adding it to its transpose, now, the upper triangle overwrites the lower.

# latentnet 2.7.1

* The "Quantile of 0" printed out for the posterior mean summary is now "2*min(Pr(>0),Pr(<0))". Also, previously, this quantity was not doubled.
* A bug in simulation for normal families have been fixed; dispersion is now handled correctly.

# latentnet 2.7.0

* A bug in the calculation of BIC has been fixed: the effective number of parameters in the likelihood was not correctly calculated. If the BIC is only used in selecting the number of clusters, the qualitative results should not be affected (i.e., the selected cluster count should be the same), and similarly if it is used to select the fixed effects to use in the model. (This BIC should not, in general, be used to select the latent space dimension or whether or not to include random effects.)
* bic.ergmm() no longer prints out a warning with random actor effects models.
* bic.ergmm() now takes an additional argument, eff.obs, specifying how to calculate the effective number of observations. The default prior to this release has been the number of actors in the network. The new default is the number of ties, as Handcock et al. (2007) recommend.
* latentnet no longer accesses unexported object in coda's namespace, as it no longer needs to redefine coda's generics.
* latentnet now exports the minimum necesary functions. Current dependencies should not be affected, but if you would like some function to be exported for your package, please let the maintainer know.
* Minor memory bugs have been fixed.

# latentnet 2.6.0

* Normal model with identity link can now estimate the dispersion (variance) parameter. (Before, it had to be specified.)
* Improved compatibility with ergm 3.2, including removal of mcmc.diagnostics() and other prototypes.
* A number of bugs have been fixed, including some code warnings from the latest versions of R.
* Kirk Li and Li Wang have been added as contributors for making many of these fixes.
* Sampson's monks dataset included with latentnet now has some additional attributes.

# latentnet 2.5.1

* A bug in importing of ergm terms for undirected networks has been fixed.

# latentnet 2.5.0

* The deprecated term `latent` has been removed. Use `euclidean` instead.
* summary() and print() functions have somewhat more intelligent defaults.
* latentnet can now use ergm() terms as covariates in the model (and the prior can still be specified). These should work seamlessly for networks without self-loops. ergmm() terms to model self-loops explicitly have been added.
* Intercepts are handled a bit more intelligently.

# latentnet 2.4.5

* A workaround for the overridden as.mcmc() generic has been added.
* All examples lines are now shorter than 100 characters.

# latentnet 2.4.4

* The implementation of Procrustes analysis has been reworked in two ways: 1) latentnet no longer depends on the shapes library; and 2) the C implementation of procrustes no longer uses the (deprecated) EISPACK, replacing it with LAPACK.
* Commonly used non-statistical functions have been extracted into a new package, statnet.common, on which latentnet depends.
* All contributors (as opposed to only the major ones) are now listed in the DESCRIPTION file.
* 3D plots can now draw segments (for undirected graphs) and arrows (for directed graphs, with package heplots) between vertices.
* More examples and tests have been made conditional on environment variable ENABLE_statnet_TESTS being set.

# latentnet 2.4.3

* In order to make upgrading from the preview release automatic, this version number has been skipped.

# latentnet 2.4.2

* Latentnet is now compliant with R 2.15.0.
* Latentnet no longer relies on a certain behavior of .C() that was changed in the development version of R.
* Many small updates, robustifications, and bug fixes have been made. In particular, crashes and weird memory behavior with recent versions of R and GCC have been fixed.
* Author and maintainer e-mail addresses and affiliations have been updated.
* License specification has been clarified.
* Latentnet now prints a warning when a matrix is passed on the LHS of the model formula instead of a network.
* Calculation of geodesic distances has been "outsourced" to the sna package, which is now a dependency.
* Fewer functions are unnecessarily masked.
* Documentation has been cleaned up and updated. Examples and tests have been shortened, speeding up R CMD check.

# latentnet 2.4.1

* Latentnet now prints a warning when computing BICs for random effects models, since we do not have a theoretical justification for that.
* A bug reported by Jorge González in the multithreaded code has been fixed.

# latentnet 2.4.0

* Removed the dependence on package 'mclust'. (The functionality used has been reimplemented internally.) This means that commercial use of the latent position clustering functionality is now possible without the UW license. The reimplementation requires 'mvtnorm' as a dependency.
* Fixed a number of bugs reported by Soriani Nicola.
* Modified the proposal tuner to propose those variables with high autocorrelation more aggressively (at the expense of the others). This should vastly improve mixing for complex models.
* Package 'coda' is now a dependency rather than a suggested package.
* Made the require()s for the suggested packages more informative.

# latentnet 2.3.0

* Added an ability to fit bilinear latent space models of Hoff (2005).
* Bug fixes, including one reported by Soriani Nicola.
* Documentation fixes.

# latentnet 2.2.3

* Fixed a bug in handling of missing data, uncovered by Xiaoyue "Maggie" Niu, and added a test for it.
* Documentation fixes.

# latentnet 2.2.2

* More documentation fixes.

# latentnet 2.2.1

* Documentation fixes.

# latentnet 2.2.0

* Ability to fit random effects (actor heterogeneity) models.
* Removed ergmm.par class --- it has been rendered obsolete by changes in R's handling of partial matching in [[ operator.
* Many bug fixes and miscellaneous changes.
