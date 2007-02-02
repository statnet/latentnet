if(!exists("mcmc.diagnostics", mode="function")){
  mcmc.diagnostics <- function(object, ...)
    UseMethod("mcmc.diagnostics")
  
  mcmc.diagnostics.default <- function(object,...){
    stop("An object must be given as an argument ")
  }
}

mcmc.diagnostics.ergmm <- function(x,burnin=FALSE,
                                   which.vars=NULL,
                                   vertex.i=c(1)){
  oldask=par("ask")
  on.exit(par(ask=oldask))
  x <- as.mcmc.list.ergmm(x,burnin,which.vars,vertex.i,vertex.i)

  autocorr.plot(x)

  par(ask=dev.interactive())
  plot(x)
  rd<-raftery.diag(x,r=0.0125)
  print(rd)
  invisible(rd)
}

as.mcmc.list.ergmm<-function(x,burnin=FALSE,
                             which.vars=NULL,
                             Z.i=c(1),
                             randeff.i=c(1)){
  n<-network.size(x$model$Yg)
  G<-x$model$G
  d<-x$model$d
  p<-x$model$p
  start<-x$control$burnin
  thin<-x$control$interval

  as.mcmc.list.ergmm.par.list(if(burnin) x$burnin else x$samples,
                              if(is.null(which.vars)) list(llk=1,
                                                           beta=1:p,
                                                           Z=cbind(rep(Z.i,each=d),rep(1:d,length(Z.i))),
                                                           sender=randeff.i,
                                                           receiver=randeff.i,
                                                           sociality=randeff.i) else which.vars,
                              start,thin)
}
