.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("latentnet",c("statnet"),FALSE)
  if(!is.null(sm)) packageStartupMessage(sm)
    
  ## Remember where we loaded this instance of latentnet, so that
  ## the snowFT slave functions could do the same.
  assign("path.to.me",file_path_as_absolute(lib),pos="package:latentnet")
  assign("nlog.double.eps",-log(.Machine[["double.eps"]]),pos="package:latentnet")

  # If the following have already been defined in the ergm package, don't duplicate. Otherwise, assign them.
  IFNOTEXISTS <- c("robust.inverse","mcmc.diagnostics","mcmc.diagnostics.default","gof","gof.default")
  for(fun in IFNOTEXISTS){
    if(!exists(fun, mode="function")){
      assign(fun, get(paste('.',fun,sep='')), pos="package:latentnet")
    }
  }
}
