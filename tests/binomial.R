#  File tests/binomial.R in package latentnet, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
library(latentnet)
data(sampson)

fit.1 <- ergmm(samplike ~ euclidean(2) + rreceiver, tofit = "mle")
eta <- function(par){
  o <- par$beta - as.matrix(dist(par$Z)) + outer(rep(1, 18), par$receiver)
  diag(o) <- -Inf
  o
}

trials <- matrix(1000 + rgeom(18*18, 1/1000), 18, 18)
diag(trials) <- 0

Ym <- matrix(rbinom(18*18, trials, plogis(eta(fit.1$mle))), 18, 18)

Y <- as.network(Ym, ignore.eval = FALSE, names.eval = "a")

wrong.trials <- matrix(floor(runif(18*18, pmin(trials, Ym), trials*2 - pmin(trials, Ym) + 1)), 18, 18)

fit.mle <- ergmm(Y ~ euclidean(2) + rreceiver, response = "a", tofit = "mle", family = "binomial", fam.par = list(trials = trials))
fit.mle.wrong <- ergmm(Y ~ euclidean(2) + rreceiver, response = "a", tofit = "mle", family = "binomial", fam.par = list(trials = wrong.trials))

stopifnot(sqrt(mean(c((eta(fit.1$mle) - eta(fit.mle$mle))^2), na.rm=TRUE)) < sqrt(mean(c((eta(fit.1$mle) - eta(fit.mle.wrong$mle))^2), na.rm=TRUE)))
stopifnot(sqrt(mean(c((plogis(eta(fit.1$mle)) - plogis(eta(fit.mle$mle)))^2), na.rm=TRUE)) < sqrt(mean(c((plogis(eta(fit.1$mle)) - plogis(eta(fit.mle.wrong$mle)))^2), na.rm=TRUE)))

fit.mkl <- ergmm(Y ~ euclidean(2) + rreceiver, response = "a", tofit = "mkl", family = "binomial", fam.par = list(trials = trials))
fit.mkl.wrong <- ergmm(Y ~ euclidean(2) + rreceiver, response = "a", tofit = "mkl", family = "binomial", fam.par = list(trials = wrong.trials))

stopifnot(sqrt(mean(c((eta(fit.1$mle) - eta(fit.mkl$mkl))^2), na.rm=TRUE)) < sqrt(mean(c((eta(fit.1$mle) - eta(fit.mkl.wrong$mkl))^2), na.rm=TRUE)))
stopifnot(sqrt(mean(c((plogis(eta(fit.1$mle)) - plogis(eta(fit.mkl$mkl)))^2), na.rm=TRUE)) < sqrt(mean(c((plogis(eta(fit.1$mle)) - plogis(eta(fit.mkl.wrong$mkl)))^2), na.rm=TRUE)))
