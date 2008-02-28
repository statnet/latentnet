library(latentnet)

data(sampson)
samp.fit <- ergmm(samplike ~ latent(d=2, G=3),control=ergmm.control(store.burnin=TRUE))

print(summary(samp.fit))
for(i in samp.fit$control$pilot.runs) mcmc.diagnostics(samp.fit,burnin=i)
mcmc.diagnostics(samp.fit)

plot(samp.fit,labels=TRUE)
plot(samp.fit,pie=TRUE)
plot(samp.fit,what="pmean")
plot(samp.fit,what="cloud")
plot(samp.fit,what="density")
plot(samp.fit,what=5)

plot(simulate(samp.fit))
plot(with(samp.fit,simulate(model,par=sample[[1]],prior=prior)))

data(tribes)
tribes.fit<-ergmm(tribes~latent(d=2,G=3),response="sign.012",family="binomial",fam.par=list(trials=2))
plot(tribes.fit,edge.col=as.matrix(tribes,"gama",m="a")*3+as.matrix(tribes,"rova",m="a")*2,pie=TRUE)

data(davis)
davis.fit<-ergmm(davis~latent(d=2,G=2))
mcmc.diagnostics(davis.fit)
plot(davis.fit,pie=TRUE)
