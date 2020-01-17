#  File R/print.ergmm.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
#' @export
print.ergmm<-function(x,...){
cat("Fitted ERGMM model:\n")
print(x$model)
}

#' @export
show.ergmm <- print.ergmm
