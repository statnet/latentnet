#  File tests/onlymcmc.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
library(statnet.common)
opttest({  # only run if ENABLE_MPI_TESTS flag is set, to avoid on windows where MPI is hard to get working. 

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

}, "parallel_MPI", testvar="ENABLE_MPI_TESTS")
