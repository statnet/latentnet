#  File R/zzz.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
# Create a place for package-wide global variables.
.latentnetEnv <- new.env()

#' @importFrom tools file_path_as_absolute
.onLoad <- function(lib, pkg){
  ## Remember where we loaded this instance of latentnet, so that
  ## the snowFT slave functions could do the same.
  .latentnetEnv$path.to.me <- file_path_as_absolute(lib)
  .latentnetEnv$nlog.double.eps <- -log(.Machine[["double.eps"]])
  .latentnetEnv$BIC.warned <- FALSE
  .RegisterKeywords()
}

#' @importFrom statnet.common statnetStartupMessage
.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("latentnet",c("statnet"),FALSE)
  if(!is.null(sm)) packageStartupMessage(sm,
                                         "NOTE: BIC calculation prior to latentnet 2.7.0 had a bug in the calculation of the effective number of parameters. See help(summary.ergmm) for details.\n",
                                         "NOTE: Prior to version 2.8.0, handling of fixed effects for directed networks had a bug: the covariate matrix was transposed.")
}

.RegisterKeywords <- function() {
  ergm_keyword(name="latent", short="latent", description="latent variable effect only used in \\pkg{latentnet}'s \\code{ergmm()}", popular=TRUE, package="latentnet")
  ergm_keyword(name="loops", short="loop", description="term is designed for networks with self-loops", popular=FALSE, package="latentnet")
}
