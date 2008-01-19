.First.lib <- function(lib, pkg){
    library.dynam("latentnet", pkg, lib)
    
    ## Remember where we loaded this instance of latentnet, so that
    ## the snowFT slave functions could do the same.
    assign("path.to.me",file_path_as_absolute(lib),pos="package:latentnet")
    assign("nlog.double.eps",-log(.Machine$double.eps),pos="package:latentnet")

    ehelp <- library(help="latentnet",lib.loc=lib,character.only=TRUE)$info[[1]]
    cat(paste(substring(ehelp[4],first=16),"\n",
              "Version ",substring(ehelp[2],first=16),
              " created on ",
              substring(ehelp[3],first=16),".\n", sep=""))
    
    cat(paste(
"copyright (c) 2003-2008, Pavel N. Krivitsky, University of Washington\n",
"                         Mark S. Handcock, University of Washington\n",
"                         and others; see LICENSE for the full list of\n",
"                         contributors.\n"
              ,sep=""))
    cat('This package is part of statnet project (http://www.statnetproject.org).\n')
    cat('For citation information type \'citation("latentnet")\'.\n')
    cat('Type \'help(package="latentnet")\' to get started.\n')
}

.Last.lib <- function(libpath){
  library.dynam.unload("latentnet",libpath)
}
