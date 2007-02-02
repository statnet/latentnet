ergmm.tuner<-function(model, start, prior, control,verbose=FALSE,start.is.list=FALSE){
  ctrl<-control
  tuning.runs<-ctrl$tuning.runs
  threads<-control$threads
  p<-model$p

  if(!start.is.list) start<-list(start)

  ctrl$samplesize<-control$tuning.runsize*with(model,
                                               p+(network.size(Yg)+1)*(d*2+sender+receiver+sociality)+1)/threads
  ctrl$interval<-1
  ctrl$burnin<-0

  next.point<-function(msjs,ldeltas){
    keep<-msjs>median(msjs)
    apply(ldeltas[keep,],2,weighted.mean,(msjs[keep]-median(msjs))+sqrt(.Machine$double.eps))
  }
    
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
  


  ldelta.start<-log(with(ctrl,c(Z.delta,Z.tr.delta,Z.scl.delta,
                                RE.delta,RE.shift.delta,
                                rep(beta.delta,length.out=p))))

  ldeltas<-matrix(numeric(0),ncol=length(ldelta.start))
  msjs<-numeric(0)
  
  for(i in 1:tuning.runs){
    if(i<floor(tuning.runs/4)) ldelta<-ldelta.start+rnorm(1,1)+rnorm(length(ldelta.start),0,1.5)
    else if(i==floor(tuning.runs/4)) ldelta<-ldelta.start
    else if(i%%4==0) ldelta<-ldelta+rnorm(1,.5)+rnorm(length(ldelta.start),0,1)
    else ldelta<-next.point(msjs,ldeltas)
    
    ldeltas<-rbind(ldeltas,ldelta)
    
    if(verbose) cat("i=",i,": delta=",paste(round(exp(ldelta),3),collapse=",")," ",sep="")
    
    msjs<-c(msjs,opt.f(ldelta))

    if(verbose) cat("msj=",msjs[length(msjs)],"\n",sep="")
  }
  

  best.delta<-exp(next.point(msjs,ldeltas))
  cat("Estimated optimal deltas=",paste(round(best.delta,3),collapse=","),"\n")
  
  list(Z.delta=best.delta[1],
       Z.tr.delta=best.delta[2],
       Z.scl.delta=best.delta[3],
       RE.delta=best.delta[4],
       RE.shift.delta=best.delta[5],
       beta.delta=best.delta[6])
}

run.proposal<-function(model, start, prior, tune.control){
  msjump(model,ergmm.MCMC.C(model,start,prior,tune.control)$samples)
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
 
msjump<-function(model,samples){

  Z.ref<-samples[[which.max(samples$llk)]]$Z
  samples<-proc.Z.mean(samples,Z.ref,FALSE)
  S<-length(samples)

  
  y<-t(sapply(1:S,function(i){
    l<-samples[[i]]

    ## Center random effects and latent space positions for the purpose of
    ## evaluating mixing.
    ## The change in density due to random effects is separated from the individual
    ## random effects. The scale of latent space positions and their centroid
    ## are also separated from individual latent space positions.

    if(model$sender){
      sender.shift<-mean(l$sender)
      l$sender<-l$sender-mean(l$sender)
    }
    if(model$receiver){
      receiver.shift<-mean(l$receiver)
      l$receiver<-l$receiver-mean(l$receiver)
    }
    if(model$sociality){
      sociality.shift<-mean(l$sociality)
      l$sociality<-l$sociality-mean(l$sociality)
    }
    if(model$d) l$Z<-scale(l$Z)
    
    o<-c(l$llk,
         pack.optim.L(l),
         if(model$sender) sender.shift,
         if(model$receiver) receiver.shift,
         if(model$sociality) sociality.shift,
         if(model$d) c(attr(l$Z,"scaled:scale"),attr(l$Z,"scaled:center")))
    o[is.nan(o)]<-0
    o
  }))

  ## Compute median square jumps.
  dy<-diff(y)
  msjs<-apply(abs(dy),2,mean)
  
  ## "Geometric Mean" criterion. Has the benefit of being insensitive
  ## to the overall scale of the parameters.
  exp(mean(log(msjs)))
}

