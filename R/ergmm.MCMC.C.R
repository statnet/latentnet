### ergmm.MCMC.C: This is a pretty minimal R function that prepares the R data to be
### passed into the C function, calls the C function to estimate the latent space model,
### and then puts the C data back into readable R storage. The hope is to separate the
### part specific to interfacing with C from the rest of the program.
### It does NOT come up with initial values --- those must be passed to it.

### If present, parameters 'sample.size' and 'interval' override those in 'control'.

### Also note that it does NOT perform the burnin. To perform the burnin,
### pass sample.size=1 and interval=burnin, and then pass the last (and only) iteration
### to the actual run.


ergmm.MCMC.C<-function(model, start, prior, control, sample.size=NULL, interval=NULL){
  Ym.noNA<- Ym <-model$Ym
  Ym.noNA[is.na(Ym.noNA)]<-0
    
  n<-network.size(model$Yg)
  d<-model$d
  G<-model$G
  p<-model$p

  if(is.null(sample.size)) sample.size<-control$sample.size
  if(is.null(interval)) interval<-control$interval
  
  if(length(prior$beta.mean)==1) prior$beta.mean<-rep(prior$beta.mean,p)
  if(length(prior$beta.var)==1) prior$beta.var<-rep(prior$beta.var,p)

  ## Figure out the design matrix.
  observed<-observed.dyads(model$Yg)
  
  if((all(observed==(diag(n)==0)) && is.directed(model$Yg)) ||
     (all(observed==lower.tri(diag(n))) && !is.directed(model$Yg)))
    observed<-NULL

  ## Sanity checks: the following block of code checks that all dimensionalities and
  ## dimensions are correct, and those optional parameters that are required by the presence
  ## of other optional parameters are present.
  
  for(i in 1:p)
    if(!all(dim(model$X[[i]])==c(n,n))) stop("Incorrect size for covariate matrices.")

  if(!is.null(start$Z)){
    if(!all(dim(start$Z)==c(n,d))) stop("Incorrect size for the starting latent positions.")
    if(is.null(control$Z.delta)) stop("Need Z-move proposal standard deviation (control$Z.delta).")
    if(G > 0){
      if(length(start$Z.K)!=n) stop("Incorrect length for the vector of starting cluster assignments.")
      if(length(start$Z.pK)!=G) stop("Incorrect length for the vector of starting cluster probabilities.")
      if(!all(dim(start$Z.mean)==c(G,d))) stop("Incorrect size for the starting cluster means.")
      if(length(start$Z.var)!=G) stop("Incorrect size for the starting cluster variances.")
    } else{
       if(length(start$Z.var)!=1) stop("Missing starting latent space variance.")
    }
  }  
  if(length(start$beta)!=p) stop("Incorrect length for the starting beta vector.")
  if(length(prior$beta.mean)!=p) stop("Incorrect length for the prior beta mean vector.")
  if(length(prior$beta.var)!=p) stop("Incorrect length for the prior beta standard deviation vector.")

  ## End Sanity checks.

  RESERVED<-2
  
#  cat("Entering C routine... ")
  Cret <- .C("ERGMM_MCMC_wrapper",
             
             sample.size=as.integer(sample.size),
             interval=as.integer(interval),
             
             n=as.integer(n),
             p=as.integer(p),
             d=as.integer(d),
             G=as.integer(G), 
             
             dir=as.integer(is.directed(model$Yg)),
             viY=as.integer(Ym.noNA),
             vdY=as.double(Ym.noNA),
             family=as.integer(model$familyID),
             iconsts=as.integer(model$iconsts),
             dconsts=as.double(model$dconsts),
             
             vX=as.double(unlist(model$X)),  
             
             llk.mcmc=double(sample.size+RESERVED),
             lpZ.mcmc=if(!is.null(start$Z))double(sample.size+RESERVED) else double(0),
             lpbeta.mcmc=if(p>0)double(sample.size+RESERVED) else double(0),
             lpLV.mcmc=if(!is.null(start$Z))double(sample.size+RESERVED) else double(0),
             
             Z=as.double(start$Z),
             
             Z.pK=if(G > 0) as.double(start$Z.pK) else double(0),
             Z.mean=if(G > 0) as.double(start$Z.mean) else double(0),
             Z.var=as.double(start$Z.var),
             Z.K=if(G > 0) as.integer(start$Z.K) else integer(0),
             
             prior.Z.var=as.double(prior$Z.var),
             prior.Z.mean.var=if(G > 0) as.double(prior$Z.mean.var) else double(0),
             prior.Z.pK=if(G > 0) as.double(prior$Z.pK) else double(0),
             prior.Z.var.df=as.double(prior$Z.var.df),
             
             Z.mcmc = double((sample.size+RESERVED)*n*d),
             Z.rate = if(d > 0) double((sample.size+RESERVED)) else double(0),
             
             K.mcmc = if(G > 0) integer(n*(sample.size+RESERVED)) else integer(0),
             Z.pK.mcmc = double(G*(sample.size+RESERVED)),
             mu.mcmc = double(d*G*(sample.size+RESERVED)),
             Z.var.mcmc = double(max(G,d>0)*(sample.size+RESERVED)),
             
             start.beta=as.double(start$beta),
             prior.beta.mean=as.double(prior$beta.mean),
             prior.beta.var=as.double(prior$beta.var),
             beta.mcmc=double((sample.size+RESERVED)*p),
             beta.rate=double((sample.size+RESERVED)),
             
             observed=as.integer(observed),

             deltas=with(control,as.numeric(c(Z.delta,group.deltas))),

             accept.all=control$accept.all,
             
             PACKAGE="latentnet")
#  cat("Finished C routine.\n")
  
  sample<-list(## MCMC Sample
                llk=Cret$llk.mcmc,
                lpZ=Cret$lpZ.mcmc,
                lpbeta=Cret$lpbeta.mcmc,
                lpLV=Cret$lpLV.mcmc,
                beta=matrix(Cret$beta.mcmc,ncol=p),
                beta.rate=Cret$beta.rate,
                Z.K = if(G>0) matrix(Cret$K.mcmc,ncol=n),
                Z.mean = if(G>0) array(Cret$mu.mcmc,dim=c((sample.size+RESERVED),G,d)),
                Z.var = if(d>0) matrix(Cret$Z.var.mcmc,ncol=max(G,1)),
                Z.pK = if(G>0) matrix(Cret$Z.pK.mcmc,ncol=G),
                Z=if(d>0)array(Cret$Z.mcmc,dim=c((sample.size+RESERVED),n,d)),
                Z.rate=if(d>0) Cret$Z.rate
                )
  class(sample)<-"ergmm.par.list"
  
  
  mcmc.mle<-sample[[1]]
  mcmc.pmode<-sample[[2]]
  sample<-del.iteration(sample,1:2)
  
  
  ## Construct the list (of lists) for return.
  out<-list(sample=sample,
            mcmc.mle=mcmc.mle,
            mcmc.pmode=mcmc.pmode)

  out
}



ergmm.MCMC.snowFT<-function(threads, reps, model.l, start.l, prior.l, control.l, sample.size.l=NULL, interval.l=NULL){
  l.sizes<-c(length(model.l),
             length(start.l),
             length(prior.l),
             length(control.l),
             length(sample.size.l),
             length(interval.l))
  l.sizes<-l.sizes[l.sizes>0]
  param.sets<-max(l.sizes)
  if(any(l.sizes!=param.sets & l.sizes!=1)) stop("Length of each input list must be either 1 or a single other number.")

  if(!require(snowFT)) stop("Package snowFT required for multithreaded MCMC!")
  mcmc.out.l<-performParallel(threads,rep(1:param.sets,reps),
                              ergmm.MCMC.snowFT.slave,
                              lib=path.to.me,
                              model.l=model.l,
                              start.l=start.l,
                              prior.l=prior.l,
                              control.l=control.l,
                              sample.size.l=sample.size.l,
                              interval.l=interval.l)
  mcmc.mle<-mcmc.out.l[[which.max(sapply(1:length(mcmc.out.l),
                                         function(i) mcmc.out.l[[i]]$mcmc.mle$llk))]]$mcmc.mle
  mcmc.pmode<-mcmc.out.l[[which.max(sapply(1:length(mcmc.out.l),
                                         function(i) mcmc.out.l[[i]]$mcmc.pmode$llk))]]$mcmc.pmode
  result.list<-list(sample=list(),mcmc.mle=mcmc.mle,mcmc.pmode=mcmc.pmode)
  for(i in 1:length(mcmc.out.l)) result.list$sample[[i]]<-mcmc.out.l[[i]]$sample
  result.list
}

ergmm.MCMC.snowFT.slave<-function(i, lib, model.l, start.l, prior.l, control.l, sample.size.l=NULL, interval.l=NULL){
  library(latentnet,lib=lib)
  ergmm.MCMC.C(model.l[[min(length(model.l),i)]],
               start.l[[min(length(start.l),i)]],
               prior.l[[min(length(prior.l),i)]],
               control.l[[min(length(control.l),i)]],
               sample.size.l[[min(length(sample.size.l),i)]],
               interval.l[[min(length(interval.l),i)]])
}

