#  File tests/nomixing.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
library(latentnet)

data(sampson)

badfit<-ergmm(samplike~euclidean(d=2,G=3)+rreceiver,control=ergmm.control(mle.maxit=3,burnin=0,interval=1,sample.size=1000,group.deltas=0,pilot.runs=0))

plot(badfit)

mcmc.diagnostics(badfit)
