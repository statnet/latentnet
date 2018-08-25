#  File R/ergmm.geodesicmatrix.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
ergmm.geodesicmatrix<-function(model){
  Ym<-model[["Ym"]]
  # For the purpose of geodesic distance, dichotomize the network about its mean.
  Ym<-Ym>mean(Ym,na.rm=TRUE)
  Ym[is.na(Ym)]<-0
  Ym <- Ym|t(Ym)
  #' @importFrom sna geodist
  geodist(Ym, count.paths=FALSE, inf.replace=nrow(Ym))$gdist
}
