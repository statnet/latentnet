### ergmm.MCMC.C: This is a pretty minimal R function that prepares the R data to be
### passed into the C function, calls the C function to estimate the latent space model,
### and then puts the C data back into readable R storage. The hope is to separate the
### part specific to interfacing with C from the rest of the program.
### It does NOT come up with initial values --- those must be passed to it.

### If present, parameters 'samplesize' and 'interval' override those in 'control'.

### Also note that it does NOT perform the burnin. To perform the burnin,
### pass samplesize=1 and interval=burnin, and then pass the last (and only) iteration
### to the actual run.


ergmm.MCMC.C<-function(model, start, prior, control, samplesize=NULL, interval=NULL){
  ## Note that passing NULL as a parameter will cause the corresponding parameter in
  ## latent_MCMC_wrapper(...) to be set to NULL when NULL is coerced to double.
  ## (as.double(NULL)==double(0))

  Ym <-model$Ym
  
  n<-network.size(model$Yg)
  d<-model$d
  G<-model$G
  p<-model$p

  if(is.null(samplesize)) samplesize<-control$samplesize
  if(is.null(interval)) interval<-control$interval
  
  if(length(prior$beta.mean)==1) prior$beta.mean<-rep(prior$beta.mean,p)
  if(length(prior$beta.var)==1) prior$beta.var<-rep(prior$beta.var,p)

  ## Figure out the design matrix.
  observed<-observed.dyads(model$Yg)

  if((observed==(diag(n)==0) && is.directed(model$Yg)) ||
     (observed==lower.tri(diag(n)) && !is.directed(model$Yg)))
    observed<-NULL

  ## Sanity checks: the following block of code checks that all dimensionalities and
  ## dimensions are correct, and those optional parameters that are required by the presence
  ## of other optional parameters are present.
  
  for(i in 1:p)
    if(!all(dim(model$X[[i]])==c(n,n))) stop("Incorrect size for covariate matrices.")

  if(!is.null(start$Z)){
    if(!all(dim(start$Z)==c(n,d))) stop("Incorrect size for the starting latent positions.")
    if(is.null(control$Z.delta)) stop("Need Z-move proposal standard deviation (control$Z.delta).")
    if(is.null(control$Z.scl.delta)) stop("Need Z-scale proposal standard deviation (control$Z.delta).")
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
             
             samplesize=as.integer(samplesize),
             interval=as.integer(interval),
             
             n=as.integer(n),
             p=as.integer(p),
             d=as.integer(d),
             G=as.integer(G), 
             
             dir=as.integer(is.directed(model$Yg)),
             viY=as.integer(Ym),
             vdY=as.double(Ym),
             family=as.integer(model$familyID),
             iconsts=as.integer(model$iconsts),
             dconsts=as.double(model$dconsts),
             
             vX=as.double(unlist(model$X)),  
             
             llk.mcmc=double(samplesize+RESERVED),
             lpZ.mcmc=if(!is.null(start$Z))double(samplesize+RESERVED) else double(0),
             lpbeta.mcmc=if(p>0)double(samplesize+RESERVED) else double(0),
             lpLV.mcmc=if(!is.null(start$Z))double(samplesize+RESERVED) else double(0),
             
             Z=as.double(start$Z),
             
             Z.pK=if(G > 0) as.double(start$Z.pK) else double(0),
             Z.mean=if(G > 0) as.double(start$Z.mean) else double(0),
             Z.var=as.double(start$Z.var),
             Z.K=if(G > 0) as.integer(start$Z.K) else integer(0),
             
             prior.Z.var=as.double(prior$Z.var),
             prior.Z.mean.var=if(G > 0) as.double(prior$Z.mean.var) else double(0),
             prior.Z.pK=if(G > 0) as.double(prior$Z.pK) else double(0),
             prior.Z.var.df=as.double(prior$Z.var.df),
             
             Z.mcmc = double((samplesize+RESERVED)*n*d),
             Z.rate = if(d > 0) double((samplesize+RESERVED)) else double(0),
             Z.rate.move.all = if(d > 0) double((samplesize+RESERVED)) else double(0),
             
             K.mcmc = if(G > 0) integer(n*(samplesize+RESERVED)) else integer(0),
             Z.pK.mcmc = double(G*(samplesize+RESERVED)),
             mu.mcmc = double(d*G*(samplesize+RESERVED)),
             Z.var.mcmc = double(max(G,d>0)*(samplesize+RESERVED)),
             
             start.beta=as.double(start$beta),
             prior.beta.mean=as.double(prior$beta.mean),
             prior.beta.var=as.double(prior$beta.var),
             beta.mcmc=double((samplesize+RESERVED)*p),
             beta.rate=double((samplesize+RESERVED)),
                          
             observed=as.integer(observed),

             deltas=with(control,as.numeric(c(Z.delta,Z.tr.delta,Z.scl.delta,
               rep(beta.delta,length.out=p)))),
             
             PACKAGE="latentnet")
#  cat("Finished C routine.\n")
  
  samples<-list(## MCMC Samples
                llk=Cret$llk.mcmc,
                lpZ=Cret$lpZ.mcmc,
                lpbeta=Cret$lpbeta.mcmc,
                lpLV=Cret$lpLV.mcmc,
                beta=matrix(Cret$beta.mcmc,ncol=p),
                beta.rate=Cret$beta.rate,
                Z.K = if(G>0) matrix(Cret$K.mcmc,ncol=n),
                Z.mean = if(G>0) array(Cret$mu.mcmc,dim=c((samplesize+RESERVED),G,d)),
                Z.var = if(d>0) matrix(Cret$Z.var.mcmc,ncol=max(G,1)),
                Z.pK = if(G>0) matrix(Cret$Z.pK.mcmc,ncol=G),
                Z=if(d>0)array(Cret$Z.mcmc,dim=c((samplesize+RESERVED),n,d)),
                Z.rate=if(d>0) Cret$Z.rate, Z.rate.move.all=if(d>0) Cret$Z.rate.move.all
                )
  class(samples)<-"ergmm.par.list"
  
  
  mcmc.mle<-samples[[1]]
  mcmc.pmode<-samples[[2]]
  samples<-del.iteration(samples,1:2)
  
  
  ## Construct the list (of lists) for return.
  out<-list(samples=samples,
            mcmc.mle=mcmc.mle,
            mcmc.pmode=mcmc.pmode)

  out
}



ergmm.MCMC.snowFT<-function(threads, reps, model.l, start.l, prior.l, control.l, samplesize.l=NULL, interval.l=NULL){
  l.sizes<-c(length(model.l),
             length(start.l),
             length(prior.l),
             length(control.l),
             length(samplesize.l),
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
                              samplesize.l=samplesize.l,
                              interval.l=interval.l)
  mcmc.mle<-mcmc.out.l[[which.max(sapply(1:length(mcmc.out.l),
                                         function(i) mcmc.out.l[[i]]$mcmc.mle$llk))]]$mcmc.mle
  
  result.list<-list(samples=list(),mcmc.mle=mcmc.mle)
  for(i in 1:length(mcmc.out.l)) result.list$samples[[i]]<-mcmc.out.l[[i]]$samples
  result.list
}

ergmm.MCMC.snowFT.slave<-function(i, lib, model.l, start.l, prior.l, control.l, samplesize.l=NULL, interval.l=NULL){
  library(latentnet,lib=lib)
  ergmm.MCMC.C(model.l[[min(length(model.l),i)]],
               start.l[[min(length(start.l),i)]],
               prior.l[[min(length(prior.l),i)]],
               control.l[[min(length(control.l),i)]],
               samplesize.l[[min(length(samplesize.l),i)]],
               interval.l[[min(length(interval.l),i)]])
}

