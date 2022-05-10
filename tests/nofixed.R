#  File tests/nofixed.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
library(latentnet)

data(sampson)

monks.nf<-ergmm(samplike~euclidean(d=2)+rreceiver-1)
mcmc.diagnostics(monks.nf)
plot(gof(monks.nf))
predict(monks.nf)
simulate(monks.nf)
print(summary(monks.nf))
