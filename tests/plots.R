library(latentnet)

data(sampson)
samp.fit <- ergmm(samplike ~ latent(d=2, G=3)+rreceiver,control=ergmm.control(store.burnin=TRUE))

print(summary(samp.fit))
for(i in samp.fit$control$pilot.runs) mcmc.diagnostics(samp.fit,burnin=i)
mcmc.diagnostics(samp.fit)

plot(samp.fit,labels=TRUE,rand.eff="receiver")
plot(samp.fit,pie=TRUE,rand.eff="receiver")
plot(samp.fit,what="pmean",rand.eff="receiver")
plot(samp.fit,what="cloud",rand.eff="receiver")
plot(samp.fit,what="density",rand.eff="receiver")
plot(samp.fit,what=5,rand.eff="receiver")

plot(simulate(samp.fit))
with(samp.fit,simulate(model,par=sample[[1]],prior=prior))
