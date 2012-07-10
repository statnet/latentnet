library(latentnet)
opttest({
rm(list=ls())
{
library(coda)
data(sampson)
monks.fit<-ergmm(samplike~euclidean(d=2,G=3))
monks.fit.mcmc<-as.mcmc.list(monks.fit)
plot(monks.fit.mcmc)
raftery.diag(monks.fit.mcmc)
}
},"as.mcmc.list.ergmm.Rd")
opttest({
rm(list=ls())
{
data(davis)
# Fit a 2D 2-cluster fit and plot.
davis.fit<-ergmm(davis~euclidean(d=2,G=2)+rsociality)
plot(davis.fit,pie=TRUE,rand.eff="sociality")
}
},"davis.Rd")
opttest({
rm(list=ls())
{
data(sampson)
## Shorter run than default.
samp.fit<-ergmm(samplike~euclidean(d=2,G=3)+rreceiver,
control=ergmm.control(burnin=1000,sample.size= 2000,interval=5))
}
},"ergmm.control.Rd")
opttest({
rm(list=ls())
{
#
# Use 'data(package = "latentnet")' to list the data sets in a
#
data(package="latentnet")
#
# Using Sampson's Monk data, lets fit a 
# simple latent position model
#
data(sampson)
samp.fit <- ergmm(samplike ~ euclidean(d=2))
#
# See if we have convergence in the MCMC
mcmc.diagnostics(samp.fit)
#
# Plot the fit
#
plot(samp.fit)
#
# Using Sampson's Monk data, lets fit a latent clustering random effects model
#
samp.fit <- ergmm(samplike ~ euclidean(d=2, G=3)+rreceiver)
#
# See if we have convergence in the MCMC
mcmc.diagnostics(samp.fit)
#
# Plot the fit.
#
plot(samp.fit)
}
},"ergmm.Rd")
opttest({
rm(list=ls())
{
#
data(sampson)
#
# test the mcmc.diagnostics function
#
gest <- ergmm(samplike ~ euclidean(d=2),
              control=ergmm.control(burnin=1000,interval=5))
summary(gest)
#
# Plot the traces and densities
#
mcmc.diagnostics(gest)
}
},"mcmc.diagnostics.ergmm.Rd")
opttest({
rm(list=ls())
{
data(sampson)
# Run two short MCMC-based fits.
samp.fit1 <- ergmm(samplike ~ euclidean(d=2, G=3),
  control=ergmm.control(burnin=1000,interval=10,sample.size=2000))
samp.fit2 <- ergmm(samplike ~ euclidean(d=2, G=3),
  control=ergmm.control(burnin=1000,interval=10,sample.size=2000))

# Combine them, and summarize the result.
samp.fit <-  merge.ergmm(samp.fit1,samp.fit2)
summary(samp.fit)
}
},"merge.ergmm.Rd")
opttest({
rm(list=ls())
{
data(sampson)
monks.fit<-ergmm(samplike~euclidean(d=2,G=3),tofit="mcmc")
heatmap(predict(monks.fit),Rowv=NA,Colv=NA)
}
},"predict.ergmm.Rd")
opttest({
rm(list=ls())
{
data(sampson)
# Fit the model for cluster sizes 1 through 4:
fits<-list(
           ergmm(samplike~euclidean(d=2,G=1)),
           ergmm(samplike~euclidean(d=2,G=2)),
           ergmm(samplike~euclidean(d=2,G=3)),
           ergmm(samplike~euclidean(d=2,G=4))
           )



# Compute the BICs for the fits and plot them:
(bics<-reshape(as.data.frame(t(sapply(fits,function(x)c(G=x$model$G,unlist(bic.ergmm(x))[c("Y","Z","overall")])))),list(c("Y","Z","overall")),idvar="G",v.names="BIC",timevar="Component",times=c("likelihood","clustering","overall"),direction="long"))

with(bics,interaction.plot(G,Component,BIC,type="b",xlab="Clusters", ylab="BIC"))

# Summarize and plot whichever fit has the lowest overall BIC:
bestG<-with(bics[bics$Component=="overall",],G[which.min(BIC)])
summary(fits[[bestG]])
plot(fits[[bestG]])
}
},"summary.ergmm.Rd")
opttest({
rm(list=ls())
{
data(tribes)
# Only model positive ties:
tribes.fit<-ergmm(tribes~euclidean(d=2,G=3),response="sign.012",family="binomial",fam.par=list(trials=2))
# Edge color must be set manually, for green ties to represent alliance
# and for red ties to represent enmity.
plot(tribes.fit,edge.col=as.matrix(tribes,"pos",m="a")*3+as.matrix(tribes,"neg",m="a")*2,pie=TRUE)
# Model both positive and negative ties:
tribes.fit3<-ergmm(tribes~euclidean(d=2,G=3),response="sign.012",family="binomial.logit",fam.par=list(trials=2))
# Edge color must be set manually, for green ties to represent alliance
# and for red ties to represent enmity.
plot(tribes.fit3,edge.col=as.matrix(tribes,"pos",m="a")*3+as.matrix(tribes,"neg",m="a")*2,pie=TRUE)
}
},"tribes.Rd")
library(latentnet)
opttest({
rm(list=ls())
{
library(coda)
data(sampson)
monks.fit<-ergmm(samplike~euclidean(d=2,G=3))
monks.fit.mcmc<-as.mcmc.list(monks.fit)
plot(monks.fit.mcmc)
raftery.diag(monks.fit.mcmc)
}
},"as.mcmc.list.ergmm.Rd")
opttest({
rm(list=ls())
{
data(davis)
# Fit a 2D 2-cluster fit and plot.
davis.fit<-ergmm(davis~euclidean(d=2,G=2)+rsociality)
plot(davis.fit,pie=TRUE,rand.eff="sociality")
}
},"davis.Rd")
opttest({
rm(list=ls())
{
data(sampson)
## Shorter run than default.
samp.fit<-ergmm(samplike~euclidean(d=2,G=3)+rreceiver,
control=ergmm.control(burnin=1000,sample.size= 2000,interval=5))
}
},"ergmm.control.Rd")
opttest({
rm(list=ls())
{
#
# Use 'data(package = "latentnet")' to list the data sets in a
#
data(package="latentnet")
#
# Using Sampson's Monk data, lets fit a 
# simple latent position model
#
data(sampson)
samp.fit <- ergmm(samplike ~ euclidean(d=2))
#
# See if we have convergence in the MCMC
mcmc.diagnostics(samp.fit)
#
# Plot the fit
#
plot(samp.fit)
#
# Using Sampson's Monk data, lets fit a latent clustering random effects model
#
samp.fit <- ergmm(samplike ~ euclidean(d=2, G=3)+rreceiver)
#
# See if we have convergence in the MCMC
mcmc.diagnostics(samp.fit)
#
# Plot the fit.
#
plot(samp.fit)
}
},"ergmm.Rd")
opttest({
rm(list=ls())
{
#
data(sampson)
#
# test the mcmc.diagnostics function
#
gest <- ergmm(samplike ~ euclidean(d=2),
              control=ergmm.control(burnin=1000,interval=5))
summary(gest)
#
# Plot the traces and densities
#
mcmc.diagnostics(gest)
}
},"mcmc.diagnostics.ergmm.Rd")
opttest({
rm(list=ls())
{
data(sampson)
# Run two short MCMC-based fits.
samp.fit1 <- ergmm(samplike ~ euclidean(d=2, G=3),
  control=ergmm.control(burnin=1000,interval=10,sample.size=2000))
samp.fit2 <- ergmm(samplike ~ euclidean(d=2, G=3),
  control=ergmm.control(burnin=1000,interval=10,sample.size=2000))

# Combine them, and summarize the result.
samp.fit <-  merge.ergmm(samp.fit1,samp.fit2)
summary(samp.fit)
}
},"merge.ergmm.Rd")
opttest({
rm(list=ls())
{
data(sampson)
monks.fit<-ergmm(samplike~euclidean(d=2,G=3),tofit="mcmc")
heatmap(predict(monks.fit),Rowv=NA,Colv=NA)
}
},"predict.ergmm.Rd")
opttest({
rm(list=ls())
{
data(sampson)
# Fit the model for cluster sizes 1 through 4:
fits<-list(
           ergmm(samplike~euclidean(d=2,G=1)),
           ergmm(samplike~euclidean(d=2,G=2)),
           ergmm(samplike~euclidean(d=2,G=3)),
           ergmm(samplike~euclidean(d=2,G=4))
           )



# Compute the BICs for the fits and plot them:
(bics<-reshape(as.data.frame(t(sapply(fits,function(x)c(G=x$model$G,unlist(bic.ergmm(x))[c("Y","Z","overall")])))),list(c("Y","Z","overall")),idvar="G",v.names="BIC",timevar="Component",times=c("likelihood","clustering","overall"),direction="long"))

with(bics,interaction.plot(G,Component,BIC,type="b",xlab="Clusters", ylab="BIC"))

# Summarize and plot whichever fit has the lowest overall BIC:
bestG<-with(bics[bics$Component=="overall",],G[which.min(BIC)])
summary(fits[[bestG]])
plot(fits[[bestG]])
}
},"summary.ergmm.Rd")
opttest({
rm(list=ls())
{
data(tribes)
# Only model positive ties:
tribes.fit<-ergmm(tribes~euclidean(d=2,G=3),response="sign.012",family="binomial",fam.par=list(trials=2))
# Edge color must be set manually, for green ties to represent alliance
# and for red ties to represent enmity.
plot(tribes.fit,edge.col=as.matrix(tribes,"pos",m="a")*3+as.matrix(tribes,"neg",m="a")*2,pie=TRUE)
# Model both positive and negative ties:
tribes.fit3<-ergmm(tribes~euclidean(d=2,G=3),response="sign.012",family="binomial.logit",fam.par=list(trials=2))
# Edge color must be set manually, for green ties to represent alliance
# and for red ties to represent enmity.
plot(tribes.fit3,edge.col=as.matrix(tribes,"pos",m="a")*3+as.matrix(tribes,"neg",m="a")*2,pie=TRUE)
}
},"tribes.Rd")
