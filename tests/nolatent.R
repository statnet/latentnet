library(latentnet)

data(sampson)

monks.nmr<-ergmm(samplike~nodematch("group")+rreceiver)
mcmc.diagnostics(monks.nmr)

print(summary(monks.nmr))
# Should produce a meaningful error message.
print(try(plot(monks.nmr)))

opttest({
  
monks.dnmr<-ergmm(samplike~nodematch("group",diff=TRUE)+rreceiver)
mcmc.diagnostics(monks.dnmr)
print(summary(monks.dnmr))

monks.dnmr2<-ergmm(samplike~nodematch("group",diff=TRUE)+rreceiver,prior=monks.dnmr$prior)
if(!all.equal(monks.dnmr2$prior,monks.dnmr$prior)) stop("Prior specification problem!")

monks.nm<-ergmm(samplike~nodematch("group"))
mcmc.diagnostics(monks.nm)
print(summary(monks.nm))

monks.dnm<-ergmm(samplike~nodematch("group",diff=TRUE))
mcmc.diagnostics(monks.dnm)
print(summary(monks.dnm))
}, "Some non-latent-space")
