### family-specific functions

# Here, nlog.double.eps=-log(.Machine[["double.eps"]]) defined in .First.lib is used to
# decide when exp(eta)==exp(eta)+1

## Bernoulli logit
lpY.Bernoulli.logit<-function(Y,eta,fam.par=NULL){
  ifelse(eta>=nlog.double.eps,eta*(Y-1),eta*Y-log1p(exp(eta)))
}
lpYc.Bernoulli.logit<-function(Y,eta,fam.par=NULL){
  ifelse(eta>=nlog.double.eps,eta*(Y-1),eta*Y-log1p(exp(eta)))
}
pY.Bernoulli.logit<-function(Y=1,eta,fam.par=NULL){
  ifelse(eta>=nlog.double.eps,exp(eta*(Y-1)),exp(eta*Y)/(exp(eta)+1))
}
dlpY.deta.Bernoulli.logit<-function(Y,eta,fam.par=NULL) Y-EY.Bernoulli.logit(eta,fam.par)
rsm.Bernoulli.logit<-function(eta,fam.par=NULL){
  n<-dim(eta)[1]
  matrix(rbinom(n*n,1,pY.Bernoulli.logit(eta=eta,fam.par=NULL)),n,n)
}
EY.Bernoulli.logit<-function(eta,fam.par=NULL) 1/(1+exp(-eta))

## Binomial logit

lpY.binomial.logit<-function(Y,eta,fam.par){
  ifelse(eta>=nlog.double.eps,eta*(Y-fam.par[["trials"]])+lchoose(fam.par[["trials"]],Y),
         eta*Y-fam.par[["trials"]]*log1p(exp(eta))+lchoose(fam.par[["trials"]],Y))
}
lpYc.binomial.logit<-function(Y,eta,fam.par){
  ifelse(eta>=nlog.double.eps,eta*(Y-fam.par[["trials"]]),
         eta*Y-fam.par[["trials"]]*log1p(exp(eta)))
}
pY.binomial.logit<-function(Y,eta,fam.par) exp(lpY.binomial.logit(Y,eta,fam.par))
dlpY.deta.binomial.logit<-function(Y,eta,fam.par) (Y-EY.binomial.logit(eta,fam.par))
rsm.binomial.logit<-function(eta,fam.par){
  n<-dim(eta)[1]
  matrix(rbinom(n*n,fam.par[["trials"]],EY.binomial.logit(eta,fam.par)/fam.par[["trials"]]),n,n)
}
EY.binomial.logit<-function(eta,fam.par) fam.par[["trials"]]/(1+exp(-eta))

## Poisson log

lpY.Poisson.log<-function(Y,eta,fam.par=NULL) dpois(Y,EY.Poisson.log(eta,fam.par),TRUE)
lpYc.Poisson.log<-function(Y,eta,fam.par=NULL) Y*eta-exp(eta)
pY.Poisson.log<-function(Y,eta,fam.par=NULL) dpois(Y,EY.Poisson.log(eta,fam.par),FALSE)
dlpY.deta.Poisson.log<-function(Y,eta,fam.par=NULL) Y-EY.Poisson.log(eta,fam.par)
rsm.Poisson.log<-function(eta,fam.par=NULL){
  n<-dim(eta)[1]
  matrix(rpois(n*n,exp(eta)),n,n)
}
EY.Poisson.log<-function(eta,fam.par=NULL) exp(eta)

## normal (linear)

lpY.normal<-function(Y,eta,fam.par=NULL) dnorm(Y,eta,sqrt(fam.par$var),TRUE)
lpYc.normal<-function(Y,eta,fam.par=NULL) -(Y-eta)^2/fam.par$var/2
pY.normal<-function(Y,eta,fam.par=NULL) dnorm(Y,eta,sqrt(fam.par$var),FALSE)
dlpY.deta.normal<-function(Y,eta,fam.par=NULL) (Y-eta)/fam.par$var
rsm.normal<-function(eta,fam.par=NULL){
  n<-dim(eta)[1]
  matrix(rnorm(n*n,eta,sqrt(fam.par$var)),n,n)
}
EY.normal<-function(eta,fam.par=NULL) eta

## Dispatcher functions

family.IDs<-list(Bernoulli=1,
                 Bernoulli.logit=1,
                 binomial=2,
                 binomial.logit=2,
                 Poisson=3,
                 Poisson.log=3,
                 Bernoulli.cont.logit=4,
                 binomial.cont.logit=5,
                 Poisson.cont.log=6,
                 normal=7)

family.names<-c("Bernoulli.logit",
                "binomial.logit",
                "Poisson.log",
                "Bernoulli.cont.logit",
                "binomial.cont.logit",
                "Poisson.cont.log",
                "normal")


lpY.fs<-c(lpY.Bernoulli.logit,
          lpY.binomial.logit,
          lpY.Poisson.log,
          lpY.Bernoulli.logit,
          lpY.binomial.logit,
          lpY.Poisson.log,
          lpY.normal)

lpYc.fs<-c(lpYc.Bernoulli.logit,
           lpYc.binomial.logit,
           lpYc.Poisson.log,
           lpYc.Bernoulli.logit,
           lpYc.binomial.logit,
           lpYc.Poisson.log,
           lpYc.normal)
          
pY.fs<-c(pY.Bernoulli.logit,
         pY.binomial.logit,
         pY.Poisson.log,
         pY.Bernoulli.logit,
         pY.binomial.logit,
         pY.Poisson.log,
         pY.normal)

dlpY.deta.fs<-c(dlpY.deta.Bernoulli.logit,
                dlpY.deta.binomial.logit,
                dlpY.deta.Poisson.log,
                dlpY.deta.Bernoulli.logit,
                dlpY.deta.binomial.logit,
                dlpY.deta.Poisson.log,
                dlpY.deta.normal)

rsm.fs<-c(rsm.Bernoulli.logit,
          rsm.binomial.logit,
          rsm.Poisson.log,
          rsm.Bernoulli.logit,
          rsm.binomial.logit,
          rsm.Poisson.log,
          rsm.normal)
          
EY.fs<-c(EY.Bernoulli.logit,
         EY.binomial.logit,
         EY.Poisson.log,
         EY.Bernoulli.logit,
         EY.binomial.logit,
         EY.Poisson.log,
         EY.normal)

fam.par.check<-function(model){
  if(model[["familyID"]]==2){
    if(is.null(model[["fam.par"]][["trials"]]))
      stop("Binomial family requires parameter `trials'.")
    model[["iconsts"]]<-model[["fam.par"]][["trials"]]
  }
  if(model[["familyID"]]==5){
    if(is.null(model[["fam.par"]][["trials"]]))
      stop("Binomial (cont) family requires parameter `trials'.")
    model[["dconsts"]]<-model[["fam.par"]][["trials"]]
  }
  if(model[["familyID"]]==7){
    if(is.null(model[["fam.par"]][["var"]]))
      stop("Binomial (cont) family requires parameter `var'.")
    model[["dconsts"]]<-model[["fam.par"]][["var"]]
  }
  model
}
