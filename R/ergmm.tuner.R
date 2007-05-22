ergmm.tuner<-function(model, start, prior, control,start.is.list=FALSE){
  ctrl<-control
  tuning.runs<-ctrl$tuning.runs
  threads<-control$threads
  p<-model$p

  if(!start.is.list) start<-list(start)

  ctrl$samplesize<-ceiling(control$tuning.runsize
                           *with(model,
                                 p+(network.size(Yg)+1)
                                 *(d*2+sender+receiver+sociality)+1)
                           /threads/ctrl$interval)
  ctrl$burnin<-0

  if(control$verbose) cat("Making",tuning.runs,"tuning runs, each with sample size",ctrl$samplesize,".\n")
  
  if(threads<=1)
    opt.f<-function(ldelta){
      my.ctrl<-ctrl
      my.ctrl[c("Z.delta","Z.tr.delta","Z.scl.delta","RE.delta","RE.shift.delta")]<-exp(ldelta[1:5])
      my.ctrl$beta.delta<-exp(ldelta[6:length(ldelta)])
      run.proposal(model,start[[1]],prior,my.ctrl)
    }
  else
    opt.f<-function(ldelta){
      my.ctrl<-ctrl
      my.ctrl[c("Z.delta","Z.tr.delta","Z.scl.delta","RE.delta","RE.shift.delta")]<-exp(ldelta[1:5])
      my.ctrl$beta.delta<-exp(ldelta[6:length(ldelta)])
      exp(mean(log(run.proposal.snowFT(threads,model,start,prior,my.ctrl))))
    }
  



  ## This is a slightly modified Complex Search (as described by Biles
  ## (1981)). I find that it works better than the gradient methods.
  
  ldelta.start<-log(with(ctrl,c(Z.delta,Z.tr.delta,Z.scl.delta,
                                RE.delta,RE.shift.delta,
                                rep(beta.delta,length.out=p))))

  gmmajs<-numeric(0)

  ldelta.min<-ldelta.start-7
  ldelta.max<-ldelta.start+7

  ldeltas<-rbind(ldelta.start,t(sapply(1:ceiling(2*length(ldelta.start)-2),function(i) runif(length(ldelta.start),ldelta.min,ldelta.max))))

  a.base<-0
  a<-0
  
  for(i in 1:tuning.runs){
    if(dim(ldeltas)[1]>length(gmmajs)) ldelta<-ldeltas[length(gmmajs)+1,]
    else{
      w<-which.min(gmmajs)
      centroid<-apply(ldeltas[-w,],2,mean)
      #a<-exp(a.base*(0.5-i/tuning.runs))/(1+exp(a.base*(0.5-i/tuning.runs))) /2 +.45
      a<-0.6
      if(i%%5==0) a<-1/a
      ldelta<-centroid+a*(centroid-ldeltas[w,])
      gmmajs<-gmmajs[-w]
      ldeltas<-rbind(ldeltas[-w,],ldelta)
    }
    
    if(control$verbose>1) cat("i=",i,": a=",round(a,1)," delta=",paste(round(exp(ldelta),3),collapse=",")," ",sep="")
    gmmaj<-opt.f(ldelta)
    gmmajs<-c(gmmajs,gmmaj)

    if(control$verbose>1) cat("gmmaj=",gmmajs[length(gmmajs)],"\n",sep="")
  }
  

  best.delta<-exp(apply(ldeltas,2,mean))
  if(control$verbose) cat("Estimated optimal deltas=",paste(round(best.delta,3),collapse=","),"\n")
  
  list(Z.delta=best.delta[1],
       Z.tr.delta=best.delta[2],
       Z.scl.delta=best.delta[3],
       RE.delta=best.delta[4],
       RE.shift.delta=best.delta[5],
       beta.delta=best.delta[6:length(best.delta)])
}

run.proposal<-function(model, start, prior, tune.control){
  gmmajump(model,ergmm.MCMC.C(model,start,prior,tune.control)$samples)
}

run.proposal.snowFT<-function(threads,model,start.l,prior,tune.control){
  unlist(performParallel(threads,1:max(threads,length(start.l)),
                         run.proposal.snowFT.slave,
                         lib=path.to.me,
                         model=model,
                         start.l=start.l,
                         prior=prior,
                         tune.control=tune.control))
}

run.proposal.snowFT.slave<-function(i,lib,model,start.l,prior,tune.control){
  library(latentnet,lib=lib)
  run.proposal(model,start.l[[min(length(start.l),i)]],prior,tune.control)
}


gmmajump<-function(model,samples){
  Z.ref<-samples[[which.max(samples$llk)]]$Z
  samples<-proc.Z.mean.C(samples,Z.ref)

  obs<-observed.dyads(model$Yg)
  y<-t(sapply(1:length(samples),function(i){
    l<-samples[[i]]

    ## Center random effects and latent space positions for the purpose of
    ## evaluating mixing.
    ## The change in density due to random effects is separated from the individual
    ## random effects. The scale of latent space positions and their centroid
    ## are also separated from individual latent space positions.

    dens<-mean(ergmm.eta.L(model,l)[obs])

    if(model$sender){
      sender.shift<-mean(l$sender)
      l$sender<-l$sender-sender.shift
    }
    if(model$receiver){
      receiver.shift<-mean(l$receiver)
      l$receiver<-l$receiver-receiver.shift
    }
    if(model$sociality){
      sociality.shift<-mean(l$sociality)
      l$sociality<-l$sociality-sociality.shift
    }
    if(model$d) l$Z<-scale(l$Z)
    
    o<-c(l$llk,
         dens,
         pack.optim.L(l),
         if(model$sender) sender.shift,
         if(model$receiver) receiver.shift,
         if(model$sociality) sociality.shift,
         if(model$d) c(attr(l$Z,"scaled:scale"),attr(l$Z,"scaled:center")))
    o[is.nan(o)]<-0
    o
  }))

  ## Compute mean squared jumps.
  dy<-diff(y)
  gmmajs<-apply(dy^2,2,trimmed.mean)
  
  
  min(gmmajs)*mean(gmmajs)

  # Since we expect (and want) a strong linear trend during the
  # burnin, we detrend before computing the acf.
#  y<-apply(y,2,function(x)lm(x~I(1:length(x)))$residuals)
  
  

  
#  ac1<-diag(acf(y,lag.max=1,plot=FALSE)$acf[2,,])

#  ac1[is.nan(ac1)]<-1

  #exp(mean(log(1-ac1)))
#  min(-ac1)
  
}
trimmed.mean<-function(x,trim=0.05){
  n<-length(x)
  mean(sort(x)[floor(trim*n):ceiling((1-trim)*n)])
}
