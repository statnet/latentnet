if(!exists("mcmc.diagnostics", mode="function")){
  mcmc.diagnostics <- function(x, ...)
    UseMethod("mcmc.diagnostics")
  
  mcmc.diagnostics.default <- function(x,...){
    stop("An object must be given as an argument ")
  }
}

mcmc.diagnostics.ergmm <- function(x,which.diags=c("cor","acf","trace","raftery"),
                                   burnin=FALSE,
                                   which.vars=NULL,
                                   vertex.i=c(1)){
  x <- as.mcmc.list.ergmm(x,burnin,which.vars,vertex.i)

  if("cor" %in% which.diags)
   novar <- apply(x,2,var)<1e-6
   if(all(novar)){
     warning("All the statistics are the same.\n")
     return(invisible())
   }else{
    colnames.x <- colnames(x)[!novar]
    x <- as.matrix(x[,!novar])
    colnames(x) <- colnames.x

    cat("\nCorrelations of sample statistics:\n")
    print(ergmm.MCMCacf(x))
  }

  if("acf" %in% which.diags)
    autocorr.plot(x)

  if("trace" %in% which.diags){
    oldask=par("ask")
    on.exit(par(ask=oldask))
    par(ask=dev.interactive())
    plot(x)
  }

  if("raftery" %in% which.diags){
    rd<-raftery.diag(x,r=0.0125)
    print(rd)
    invisible(rd)
  }
}

as.mcmc.list.ergmm<-as.mcmc.ergmm<-function(x,burnin=FALSE,
                             which.vars=NULL,
                             vertex.i=c(1)){
  n<-network.size(x$model$Yg)
  G<-x$model$G
  d<-x$model$d
  p<-x$model$p
  start<-x$control$burnin
  thin<-x$control$interval

  as.mcmc.list.ergmm.par.list(if(burnin) x$burnin else x$samples,
                              if(is.null(which.vars)) list(llk=1,
                                                           beta=1:p,
                                                           Z=cbind(rep(vertex.i,each=d),rep(1:d,length(vertex.i))),
                                                           sender=vertex.i,
                                                           receiver=vertex.i,
                                                           sociality=vertex.i) else which.vars,
                              start,thin)
}
