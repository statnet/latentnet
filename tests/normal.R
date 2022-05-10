#  File tests/normal.R in package latentnet, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
library(statnet.common)
opttest({
library(latentnet)
n<-20

y.var<-1/16
y.var.prior<-4*y.var

y<-as.network(matrix(1,n,n),dir=TRUE)

cat("Euclidean:\n")

Z<-rbind(cbind(rnorm(n/2,3),rnorm(n/2,3)),
         cbind(rnorm(n/2,-3),rnorm(n/2,-3)))

dm<--as.matrix(dist(Z))

ym<-rnorm(n*n,2+dm,sqrt(y.var))

set.edge.value(y,"v",ym)
image(as.matrix(y,a="v",m="a"))
y.fit<-ergmm(y~euclidean(d=2,G=2),response="v",family="normal",fam.par=list(prior.var=y.var.prior,prior.var.df=2),verbose=TRUE)

Z.mkl<-plot(y.fit,Z.ref=Z)
points(Z,pch=5)
cat("Mean squared difference:",sum((Z.mkl-Z)^2),"\n")


# Simulate from the fit.
y.sim <- simulate(y.fit, nsim=1)

cat("Bilinear:\n")

dm<-tcrossprod(Z)

ym<-rnorm(n*n,2+dm,sqrt(y.var))

set.edge.value(y,"v",ym)
image(as.matrix(y,a="v",m="a"))
y.fit<-ergmm(y~bilinear(d=2,G=2),response="v",family="normal",fam.par=list(prior.var=y.var.prior,prior.var.df=2),verbose=TRUE)

Z.mkl<-plot(y.fit,Z.ref=Z)
points(Z,pch=5)
cat("Mean squared difference:",sum((Z.mkl-Z)^2),"\n")


cat("Squared Euclidean:\n")

dm<--as.matrix(dist(Z))^2

ym<-rnorm(n*n,2+dm,sqrt(y.var))

set.edge.value(y,"v",ym)
image(as.matrix(y,a="v",m="a"))
y.fit<-ergmm(y~euclidean2(d=2,G=2),response="v",family="normal",fam.par=list(prior.var=y.var.prior,prior.var.df=2),verbose=TRUE)

Z.mkl<-plot(y.fit,Z.ref=Z)
points(Z,pch=5)
cat("Mean squared difference:",sum((Z.mkl-Z)^2),"\n")


# Simulate from the fit.
y.sim <- simulate(y.fit, nsim=1)

cat("No latent space:\n")

ym<-rnorm(n*n,0,sqrt(y.var))
set.edge.value(y,"v",ym)
image(as.matrix(y,a="v",m="a"))
y.fit<-ergmm(y~1,response="v",family="normal",fam.par=list(prior.var=y.var.prior,prior.var.df=2),verbose=TRUE)
summary(y.fit)


},"Normal response variable")
