get.group.deltas<-function(group.deltas, model, pilot.samples, control){
  if(!is.matrix(control$group.deltas)){
    if(!is.null(pilot.samples)){
      beta.ext<-cbind(if(model$p) pilot.samples$beta, # covariate coefs
                      if(model$d) log(apply(sqrt(apply(pilot.samples$Z^2,1:2,sum)),1,mean)), # scale of Z
                      if(model$d) apply(pilot.samples$Z,c(1,3),mean), # center of Z
                      if(model$sender) apply(pilot.samples$sender,1,mean), # sender eff.
                      if(model$receiver) apply(pilot.samples$receiver,1,mean), # receiver eff.
                      if(model$sociality) apply(pilot.samples$sociality,1,mean)) # sociality eff.
      cov.beta.ext<-cov(beta.ext)*control$group.deltas.mul


      ## Zero-out those rows in the variance-covariance matrix that are to be proposed independently.
      cov.pos<-0
      if(model$p){
        if("beta" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+(1:model$p))
        cov.pos<-cov.pos+model$p
      }
      
      if(model$d){
        if("Z.scl" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+1)
        if("Z.tr" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+1+(1:model$d))
        cov.pos<-cov.pos+1+model$d
      }

      if(model$sender){
        if("sender" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+1)
        cov.pos<-cov.pos+1
      }

      if(model$receiver){
        if("receiver" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+1)
        cov.pos<-cov.pos+1
      }
      
      if(model$sociality){
        if("sociality" %in% control$propose.ind) cov.beta.ext<-partial.diag(cov.beta.ext,cov.pos+1)
        cov.pos<-cov.pos+1
      }

      ## Take the Choletsky decomposition of the empirical covariance matrix.
      control$group.deltas<-chol(cov.beta.ext)
      
    }else if(length(control$group.deltas)==1){
      nterms<-model$p+model$d+(if(model$d) 1 else 0)+
        model$sender+model$receiver+model$sociality
      control$group.deltas<-1/sapply(1:model$p,function(i) sqrt(mean((model$X[[i]][observed.dyads(model$Yg)])^2)))
      if(model$d) control$group.deltas<-c(control$group.deltas, 0.05,
                                          rep(control$Z.delta,model$d))
      control$group.deltas<-c(control$group.deltas,rep(1,model$sender+model$receiver+model$sociality))
      control$group.deltas<-diag(control$group.deltas*control$group.deltas.mul*2/(1+nterms),nrow=nterms)
      
    }else if(length(control$group.deltas)>1){
      control$group.deltas<-diag(control$group.deltas)
    }
  }

  control
}

## Zero the off-diagonal entries in a particular set of rows and columns.
partial.diag<-function(m,idx){
  m<-as.matrix(m)
  saved<-diag(m)
  m[idx,]<-0
  m[,idx]<-0
  diag(m)<-saved
  m
}
