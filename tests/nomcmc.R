#  File tests/nomcmc.R in package latentnet, part of the
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

mleonly<-ergmm(samplike~euclidean(d=2),tofit="mle")

# Should skip MCMC.
if(!is.null(mleonly$sample)) stop("MCMC should not be run for MLE!")

# Should give an informative error message.
print(try(plot(mleonly)))

# Should plot OK.
plot(mleonly,what="mle")

# Should give an informative error message.
print(try(mcmc.diagnostics(mleonly)))
