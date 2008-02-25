library(latentnet)

data(sampson)

badfit<-ergmm(samplike~latent(d=2,G=3)+rreceiver,control=ergmm.control(mle.maxit=3,burnin=0,interval=1,group.deltas=0,pilot.runs=0))

plot(badfit)

mcmc.diagnostics(badfit)
