.onAttach <- function(lib, pkg){
  info <- packageDescription("latentnet")
    
  ## Remember where we loaded this instance of latentnet, so that
  ## the snowFT slave functions could do the same.
  assign("path.to.me",file_path_as_absolute(lib),pos="package:latentnet")
  assign("nlog.double.eps",-log(.Machine[["double.eps"]]),pos="package:latentnet")

  # If the following have already been defined in the ergm package, don't duplicate. Otherwise, assign them.
  IFNOTEXISTS <- c("robust.inverse","mcmc.diagnostics","mcmc.diagnostics.default","gof")
  for(fun in IFNOTEXISTS){
    if(!exists(fun, mode="function")){
      assign(fun, get(paste('.',fun,sep='')), pos="package:latentnet")
    }
  }

  packageStartupMessage(paste(
                          "copyright (c) 2003-2009, Pavel N. Krivitsky, Pennsylvania State University\n",
                          "                         Mark S. Handcock, University of California, Los Angeles\n",
                          "                         and others; see LICENSE for the full list of\n",
                          "                         contributors.\n",
                          'Based on "statnet" project software (statnet.org).\n',
                          'For license and citation information see statnet.org/attribution\n',
                          'or type citation("ergm").\n', sep="")
                        )
}
