######################################################################
#
# zzz.r
#
# copyright (c) 2003-2007, Pavel N. Krivitsky, University of Washington
#                          Mark S. Handcock, University of Washington
#                          David R. Hunter, Penn State University
#                          Carter T. Butts, University of California - Irvine
#                          Martina Morris, University of Washington
# written December 2003 edited afterwards
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the latentnet statnet package
#
# .First.lib is run when the package is loaded with library(latentnet)
#
######################################################################

.First.lib <- function(lib, pkg){
    library.dynam("latentnet", pkg, lib)
    
    ## Remember where we loaded this instance of latentnet, so that
    ## the snowFT slave functions could do the same.
    assign("path.to.me",file_path_as_absolute(lib),pos="package:latentnet")
    assign("nlog.double.eps",-log(.Machine$double.eps),pos="package:latentnet")
    
    cat(paste(
"copyright (c) 2003-2007, Pavel N. Krivitsky, University of Washington\n",
"                         Mark S. Handcock, University of Washington\n",
"                         Susan Shortreed, University of Washington\n",
"                         Jeremy Tantrum, University of Washington\n",
"                         Peter Hoff, University of Washington\n",sep=""))
    cat('See http://www.csde.washington.edu/statnet\n')
    cat('Type help(package="latentnet") to get started.\n')
}

.Last.lib <- function(libpath){
  library.dynam.unload("latentnet",libpath)
}
