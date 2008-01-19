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
                                   vertex.i=c(1),...){
  extraneous.argcheck(...)
  x <- as.mcmc.list.ergmm(x,burnin,which.vars,vertex.i)
  oldask=par("ask")
  on.exit(par(ask=oldask))
  par(ask=dev.interactive())

  if("cor" %in% which.diags){
    x.ac<-autocorr(x,lag=0:1)
    for(chain in seq(along=x.ac)){
      cat(paste("Chain",chain,"\n"))
      for(i in 1:2){
        cat(paste(dimnames(x.ac[[chain]])[[1]][i],"\n"))
        print(x.ac[[chain]][i,,])
        cat("\n")
      }
    }
  }

  if("acf" %in% which.diags)
    autocorr.plot(x)

  if("trace" %in% which.diags){
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
                                                           Z=cbind(rep(vertex.i,each=d),rep(1:d,length(vertex.i))))
                              else which.vars,
                              start,thin)
}
