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
plot(with(samp.fit,simulate(model,par=sample[[1]],prior=prior)))

data(tribes)
tribes.fit<-ergmm(tribes~latent(d=2,G=3),response="sign.012",family="binomial",fam.par=list(trials=2))
plot(tribes.fit,edge.col=as.matrix(tribes,"pos",m="a")*3+as.matrix(tribes,"neg",m="a")*2,pie=TRUE)

data(davis)
davis.fit<-ergmm(davis~latent(d=2,G=2)+rsociality)
mcmc.diagnostics(davis.fit)
plot(davis.fit,pie=TRUE,rand.eff="sociality")

tribes.fit3<-ergmm(tribes~latent(d=3,G=3),response="sign.012",family="binomial",fam.par=list(trials=2))
plot(tribes.fit3,pie=TRUE)
