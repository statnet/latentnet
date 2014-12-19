#  File R/zzz.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("latentnet",c("statnet"),FALSE)
  if(!is.null(sm)) packageStartupMessage(sm)
    
  ## Remember where we loaded this instance of latentnet, so that
  ## the snowFT slave functions could do the same.
  assign("path.to.me",file_path_as_absolute(lib),pos="package:latentnet")
  assign("nlog.double.eps",-log(.Machine[["double.eps"]]),pos="package:latentnet")
}
