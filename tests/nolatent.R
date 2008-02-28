library(latentnet)

data(sampson)

monks.nm<-ergmm(samplike~nodematch("group"))
mcmc.diagnostics(monks.nm)

print(summary(monks.nm))
# Should produce a meaningful error message.
print(try(plot(monks.nm)))

monks.dnm<-ergmm(samplike~nodematch("group",diff=TRUE))
mcmc.diagnostics(monks.dnm)
print(summary(monks.dnm))

monks.dnm2<-ergmm(samplike~nodematch("group",diff=TRUE),prior=monks.dnm$prior)
if(!all.equal(monks.dnm2$prior,monks.dnm$prior)) stop("Prior specification problem!")
