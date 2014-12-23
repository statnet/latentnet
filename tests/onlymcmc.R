#  File tests/onlymcmc.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
library(latentnet)

# Also test parallel code.
data(sampson)
onlymcmc <- ergmm(samplike ~ euclidean(d=2, G=3)+rreceiver,tofit="mcmc",control=control.ergmm(burnin=100,interval=1,sample.size=1000,pilot.runs=0,threads=2))

mcmc.diagnostics(onlymcmc)

# Should give an informative error message.
print(try(plot(onlymcmc)))

plot(simulate(onlymcmc))
plot(with(onlymcmc,simulate(model,par=sample[[1]],prior=prior)))

heatmap(predict(onlymcmc),Rowv=NA,Colv=NA)
