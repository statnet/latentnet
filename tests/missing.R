#  File tests/missing.R in package latentnet, part of the
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
ym<-as.matrix(samplike)
ym[1,5]<-NA
yg<-as.network(ym,matrix.type="adjacency")
samp.fit<-ergmm(yg~euclidean(d=2))

print(summary(samp.fit))
plot(samp.fit,labels=TRUE)
