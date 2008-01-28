get.group.deltas<-function(group.deltas, model, pilot.samples, control){
  if(!is.matrix(control$group.deltas)){
    if(!is.null(pilot.samples)){
      beta.ext<-cbind(if(model$p) pilot.samples$beta, # covariate coefs
##                      if(model$d) log(apply(sqrt(apply(pilot.samples$Z^2,1:2,sum)),1,mean)), # scale of Z
                      if(model$d) apply(pilot.samples$Z,c(1,3),mean), # center of Z
                      if(model$sender) apply(pilot.samples$sender,1,mean), # sender eff.
                      if(model$receiver) apply(pilot.samples$receiver,1,mean), # receiver eff.
                      if(model$sociality) apply(pilot.samples$sociality,1,mean)) # sociality eff.
      control$group.deltas<-chol(cov(beta.ext)*control$group.deltas.mul)
    }else if(length(control$group.deltas)==1){
      nterms<-model$p+model$d+#(if(model$d) 1 else 0)+
        model$sender+model$receiver+model$sociality
      control$group.deltas<-1/sapply(1:model$p,function(i) sqrt(mean((model$X[[i]][observed.dyads(model$Yg)])^2)))
      if(model$d) control$group.deltas<-c(control$group.deltas,## 0.05,
                                          rep(control$Z.delta,model$d))
      control$group.deltas<-c(control$group.deltas,rep(1,model$sender+model$receiver+model$sociality))
      control$group.deltas<-diag(control$group.deltas*control$group.deltas.mul*2/(1+nterms),nrow=nterms)
    }else if(length(control$group.deltas)>1){
      control$group.deltas<-diag(control$group.deltas)
    }
  }

  control
}
