### family-specific functions

# Here, nlog.double.eps=-log(.Machine$double.eps) defined in .First.lib is used to
# decide when exp(eta)==exp(eta)+1

## Bernoulli logit
lpY.Bernoulli.logit<-function(Y,eta,fam.par=NULL){
  ifelse(eta>=nlog.double.eps,eta*(Y-1),eta*Y-log(exp(eta)+1))
}
lpYc.Bernoulli.logit<-function(Y,eta,fam.par=NULL){
  ifelse(eta>=nlog.double.eps,eta*(Y-1),eta*Y-log(exp(eta)+1))
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
  ifelse(eta>=nlog.double.eps,eta*(Y-fam.par$trials)+lchoose(fam.par$trials,Y),
         eta*Y-fam.par$trials*log(exp(eta)+1)+lchoose(fam.par$trials,Y))
}
lpYc.binomial.logit<-function(Y,eta,fam.par){
  ifelse(eta>=nlog.double.eps,eta*(Y-fam.par$trials),
         eta*Y-fam.par$trials*log(exp(eta)+1))
}
pY.binomial.logit<-function(Y,eta,fam.par) exp(lpY.binomial.logit(Y,eta,fam.par))
dlpY.deta.binomial.logit<-function(Y,eta,fam.par) (Y-EY.binomial.logit(eta,fam.par))
rsm.binomial.logit<-function(eta,fam.par){
  n<-dim(eta)[1]
  matrix(rbinom(n*n,fam.par$trials,EY.binomial.logit(eta,fam.par)/fam.par$trials),n,n)
}
EY.binomial.logit<-function(eta,fam.par) fam.par$trials/(1+exp(-eta))

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

mk.lpY.f<-function(family="Bernoulli"){
  switch(family,
         Bernoulli=lpY.Bernoulli.logit,
         binomial=lpY.binomial.logit,
         Poisson=lpY.Poisson.log,
         stop(paste("Family",family,"is not supported.")))
}

mk.lpYc.f<-function(family="Bernoulli"){
  switch(family,
         Bernoulli=lpYc.Bernoulli.logit,
         binomial=lpYc.binomial.logit,
         Poisson=lpYc.Poisson.log,
         stop(paste("Family",family,"is not supported.")))
}

mk.dlpY.deta.f<-function(family="Bernoulli"){
  switch(family,
         Bernoulli=dlpY.deta.Bernoulli.logit,
         binomial=dlpY.deta.binomial.logit,
         Poisson=dlpY.deta.Poisson.log,
         stop(paste("Family",family,"is not supported.")))
}

mk.pY.f<-function(family="Bernoulli"){
  switch(family,
         Bernoulli=pY.Bernoulli.logit,
         binomial=pY.binomial.logit,
         Poisson=pY.Poisson.log,
         stop(paste("Family",family,"is not supported.")))
}

mk.EY.f<-function(family="Bernoulli"){
  switch(family,
         Bernoulli=EY.Bernoulli.logit,
         binomial=EY.binomial.logit,
         Poisson=EY.Poisson.log,
         stop(paste("Family",family,"is not supported.")))
}

mk.rsm.f<-function(family="Bernoulli"){
  switch(family,
         Bernoulli=rsm.Bernoulli.logit,
         binomial=rsm.binomial.logit,
         Poisson=rsm.Poisson.log,
         stop(paste("Family",family,"is not supported.")))
}


getYm<-function(Yg,response=NULL){
  if(is.null(response)){
    return(as.matrix.network(Yg, matrix.type="adjacency"))
  }else{
    if(is.matrix(response)) return(response)
    else return(as.matrix.network(Yg, response, matrix.type="adjacency"))
  }
}
  

ergmm.eta<-function(Yg,beta=NULL,X=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL){
  n<-Yg$gal$n
  dir<-Yg$gal$dir

  eta<-matrix(0,n,n)
  
  if(!is.null(Z))
    eta<-eta-as.matrix(dist(Z))

  if(!is.null(beta))
    for(k in 1:length(beta))
      eta<-eta+beta[k]*X[[k]]

  if(!is.null(sociality)) sender=receiver=sociality
  
  if(!is.null(sender)){
    eta<-eta+sender
  }

  if(!is.null(receiver))
    eta<-t(t(eta)+receiver)

  return(eta)
}

ergmm.eta.L<-function(model,theta){
  ergmm.eta(model$Yg,
            theta$beta,
            model$X,
            theta$Z,
            theta$sender,
            theta$receiver,
            theta$sociality)
}

ergmm.EY.L<-function(model,theta){
  eta<-ergmm.eta.L(model,theta)
  mk.EY.f(model$family)(eta)
}

ergmm.loglike<-function(Yg,response=NULL,Ym=NULL,family="Bernoulli",fam.par=NULL,beta=NULL,X=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL,partial.eta=NULL,up.to.const=FALSE,flyapart.penalty=0){
  if(is.null(partial.eta))
    eta<-matrix(0,Yg$gal$n,Yg$gal$n)
  else eta<-partial.eta

  Ym <- if(is.null(Ym))getYm(Yg,response) else Ym
  
  eta<-eta+ergmm.eta(Yg,beta,X,Z,sender,receiver,sociality)
  obs<-observed.dyads(Yg)
  lpY<-if(up.to.const) mk.lpYc.f(family)(Ym[obs],eta[obs],fam.par) else mk.lpY.f(family)(Ym[obs],eta[obs],fam.par)
  if(flyapart.penalty)
    return(sum(lpY)-flyapart.penalty*(sum(Z^2)+
                                      sum(beta^2)+
                                      sum(sender^2)+
                                      sum(receiver^2)+
                                      sum(sociality^2)))
  return(sum(lpY))
}

ergmm.loglike.grad<-function(Yg,response=NULL,Ym=NULL,family="Bernoulli",fam.par=NULL,beta=NULL,X=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL,partial.eta=NULL,flyapart.penalty=0){
  if(is.null(partial.eta))
    eta<-matrix(0,Yg$gal$n,Yg$gal$n)
  else eta<-partial.eta
  Ym <- if(is.null(Ym))getYm(Yg,response) else Ym
  obs<-observed.dyads(Yg)
  eta<-eta+ergmm.eta(Yg,beta,X,Z,sender,receiver,sociality)
  n<-network.size(Yg)

  dlpY.deta <- mk.dlpY.deta.f(family)(Ym,eta,fam.par)

  grad<-list()
  
  if(!is.null(beta)) grad$beta <- sapply(1:length(beta),function(k) sum(dlpY.deta*X[[k]]*obs))

  if(!is.null(Z)){
    d<-dim(Z)[2]
    Z.invdist<- -as.matrix(dist(Z))
    Z.invdist[Z.invdist==0]<-Inf
    Z.invdist<-1/Z.invdist

    grad$Z<-matrix(0,n,d)
    for(k in 1:d){
      Z.k.d<-sapply(1:n,function(i) sapply(1:n, function(j) Z[j,k]-Z[i,k]))*Z.invdist
      for(i in 1:n)
        for(j in 1:n)
          grad$Z[i,k]<-grad$Z[i,k]+(Z.k.d[i,j]*dlpY.deta[i,j]*obs[i,j]-
                                    Z.k.d[j,i]*dlpY.deta[j,i]*obs[j,i])
    }
  }

  if(!is.null(sociality))
    grad$sociality <- sapply(1:n,function(i) sum(dlpY.deta[i,]*obs[i,])+sum(dlpY.deta[,i]*obs[,i]))
  else{
    if(!is.null(sender)) grad$sender <- sapply(1:n,function(i) sum(dlpY.deta[i,]*obs[i,]))
    if(!is.null(receiver)) grad$receiver <-  sapply(1:n,function(i) sum(dlpY.deta[,i]*obs[,i]))
  }

  if(flyapart.penalty){
    if(!is.null(grad$Z)) grad$Z<-grad$Z-flyapart.penalty*2*Z
    if(!is.null(grad$beta)) grad$beta<-grad$beta-flyapart.penalty*2*beta
    if(!is.null(grad$sender)) grad$sender<-grad$sender-flyapart.penalty*2*sender
    if(!is.null(grad$receiver)) grad$receiver<-grad$receiver-flyapart.penalty*2*receiver
    if(!is.null(grad$sociality)) grad$sociality<-grad$sociality-flyapart.penalty*2*sociality
  }
  
  grad
}
  
  

ergmm.loglike.C<-function(Yg,response=NULL,Ym=NULL,family="Bernoulli",iconsts=NULL,dconsts=NULL,beta=NULL,X=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL){
  Y <-if(is.null(Ym)) getYm(Yg,response) else Ym
  n<-Yg$gal$n
  k<-if(!is.null(Z)) dim(Z)[2] else 0
  p<-length(X)

  ## Figure out the design matrix.
  observed<-observed.dyads(Yg)

  if((observed==(diag(n)==0) && is.directed(Yg)) ||
     (observed==lower.tri(diag(n)) && !is.directed(Yg)))
    observed<-NULL

  familyID<-switch(family,
                 Bernoulli=0,
                 binomial=1,
                 Poisson=2)

  ## Sanity checks: the following block of code checks that all dimensionalities and
  ## dimensions are correct, and those optional parameters that are required by the presence
  ## of other optional parameters are present.

  if(familyID==1 && length(iconsts)!=1)
    stop("Binomial family requires parameter n.")
  
  for(i in 1:p)
    if(!all(dim(X[[i]])==c(n,n))) stop("Incorrect size for covariate matrices.")

  if(!is.null(Z)){
    if(!all(dim(Z)==c(n,k))) stop("Incorrect size for the latent positions.")
  }  
  if(length(beta)!=p) stop("Incorrect length for the beta vector.")

  if(!is.null(sociality)){
    if(length(sender)!=n) stop("Incorrect length for the vector of sociality effects.")
  }
  if(!is.null(sender)){
    if(length(sender)!=n) stop("Incorrect length for the vector of sender effects.")
  }
  if(!is.null(receiver)){
    if(length(receiver)!=n) stop("Incorrect length for the vector of receiver effects.")
  }
  ## End Sanity checks.
  
  ret <- .C("ERGMM_lp_Y_wrapper",
            n=as.integer(n), p=as.integer(p),
            d=as.integer(k),
            
            dir=as.integer(is.directed(Yg)),
            vY=as.integer(Y),
            family=as.integer(familyID),iconsts=as.integer(iconsts),dconsts=as.integer(dconsts),
            
            vX=as.double(unlist(X)),  
            
            Z=as.double(Z),
            
            beta=as.double(beta),

            sender=if(is.null(sociality))as.double(sender) else as.double(sociality),
            receiver=as.double(receiver), lock.RE=!is.null(sociality),
            
            observed=as.integer(observed),

            llk=double(1),
            PACKAGE="latentnet")
  

  ret$llk
}

observed.dyads<-function(Yg){
  observed.dyads<-get.network.attribute(Yg,"design")
  if(is.null(observed.dyads))
    observed.dyads<-matrix(TRUE,Yg$gal$n,Yg$gal$n)
  else
    observed.dyads<-as.matrix.network(observed.dyads,matrix.type="adjacency")==0
      
  if(!is.directed(Yg)) observed.dyads[upper.tri(observed.dyads)]<-FALSE

  if(!has.loops(Yg)) diag(observed.dyads)<-FALSE
  
  observed.dyads
}

pack.optim<-function(beta=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL,
                     Z.var=NULL,Z.mean=NULL,
                     sender.var=NULL,receiver.var=NULL,sociality.var=NULL,fit.vars=NULL){
  if(is.null(fit.vars))
    return(c(beta,Z,sender,receiver,sociality,
             Z.var,Z.mean,
             sender.var,receiver.var,sociality.var))
  else
    return(c(if(fit.vars$beta)beta,if(fit.vars$Z)Z,
             if(fit.vars$sender)sender,if(fit.vars$receiver)receiver,if(fit.vars$sociality)sociality,
             if(fit.vars$Z.var)Z.var,if(fit.vars$Z.mean)Z.mean,
             if(fit.vars$sender.var)sender.var,if(fit.vars$receiver.var)receiver.var,if(fit.vars$sociality.var)sociality.var))
}
pack.optim.L<-function(theta,fit.vars=NULL){
  c(theta$beta,theta$Z,theta$sender,theta$receiver,theta$sociality,
    theta$Z.var,theta$Z.mean,
    theta$sender.var,theta$receiver.var,theta$sociality.var,fit.vars=fit.vars)
}

reg.fit.vars<-function(fit.vars){
  for(name in c("beta","Z","sender","receiver","sociality",
                "Z.var","Z.mean",
                "sender.var","receiver.var","sociality.var"))
    if(!(name %in% names(fit.vars)))fit.vars[[name]]<-FALSE
  fit.vars
}

inv.fit.vars<-function(fit.vars){
  for(name in c("beta","Z","sender","receiver","sociality",
                "Z.var","Z.mean",
                "sender.var","receiver.var","sociality.var"))
    if(name %in% names(fit.vars))
      fit.vars[[name]]<-!fit.vars[[name]]
  fit.vars
}

FIT_ALL<-list(beta=TRUE,Z=TRUE,sender=TRUE,receiver=TRUE,sociality=TRUE,
              Z.var=TRUE,Z.mean=TRUE,
              sender.var=TRUE,receiver.var=TRUE,sociality.var=TRUE)

FIT_MLE<-list(beta=TRUE,Z=TRUE,sender=TRUE,receiver=TRUE,sociality=TRUE)

unpack.optim<-function(v,fit.vars,n,p=NULL,d=NULL,G=NULL){
  if(is.null(p)) p<-0
  if(is.null(d)) d<-0
  if(is.null(G)) G<-0
  attach(fit.vars)
  v.must.be<-(beta*p +
              Z*n*d +
              sender*n +
              receiver*n +
              sociality*n +
              Z.var*(d>0)*max(1,G) +
              Z.mean*G*d +
              sender.var +
              receiver.var +
              sociality.var)
  if(length(v)!=v.must.be){
    detach(fit.vars)
    stop(paste("Input vector wrong length: ", length(v),
               " but should be ",v.must.be,".",
               sep=""))
  }
  pos<-0
  ret<-list()
  if(beta && p>0){
    ret$beta<-v[pos+1:p]
    pos<-pos+p
  }

  if(Z && d>0){
    ret$Z<-matrix(v[pos+1:(n*d)],nrow=n,ncol=d)
    pos<-pos+n*d
  }

  if(sender){
    ret$sender<-v[pos+1:n]
    pos<-pos+n
  }

  if(receiver){
    ret$receiver<-v[pos+1:n]
    pos<-pos+n
  }

  if(sociality){
    ret$sociality<-v[pos+1:n]
    pos<-pos+n
  }

  if(Z.var && d>0){
    ret$Z.var<-v[pos+1:max(1,G)]
    pos<-pos+max(1,G)
  }
  
  if(Z.mean && d>0 && G>0){
    ret$Z.mean<-matrix(v[pos+1:(G*d)],nrow=G,ncol=d)
    pos<-pos+G*d
  }
  
  if(sender.var){
    ret$sender.var<-v[pos+1]
    pos<-pos+1
  }
  
  if(receiver.var){
    ret$receiver.var<-v[pos+1]
    pos<-pos+1
  }

  if(sociality.var){
    ret$sociality.var<-v[pos+1]
    pos<-pos+1
  }
  detach(fit.vars)
  
  ret
}

mk.lp.optim.fs<-function(fit.vars,
                         Yg,response=NULL,Ym=NULL,family="Bernoulli",fam.par=NULL,
                         beta=NULL,X=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL,
                         Z.var=NULL,Z.mean=NULL,Z.K=NULL,
                         sender.var=NULL,receiver.var=NULL,sociality.var=NULL,
                         prior=NULL,opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV"),
                         partial.eta=NULL,flyapart.penalty=0){
  n<-network.size(Yg)
  p<-length(beta)
  d<- if(is.null(Z)) 0 else dim(Z)[2]
  G<-if(is.null(Z.mean)) 0 else dim(Z.mean)[1]
  fit.vars<-reg.fit.vars(fit.vars)

  return(list(
              f=function(v){
                V<-unpack.optim(v,fit.vars,n,
                                p=p,
                                d=d,
                                G=G)
                ergmm.lp(Yg,response,Ym=Ym,family=family,fam.par=fam.par,
                         beta=if(fit.vars$beta) V$beta else beta,
                         X=X,
                         Z=if(fit.vars$Z) V$Z else Z,
                         sender=if(fit.vars$sender) V$sender else sender,
                         receiver=if(fit.vars$receiver) V$receiver else receiver,
                         sociality=if(fit.vars$sociality) V$sociality else sociality,
                         Z.var=if(fit.vars$Z.var) V$Z.var else Z.var,
                         Z.mean=if(fit.vars$Z.mean) V$Z.mean else Z.mean,
                         Z.K=Z.K,
                         sender.var=if(fit.vars$sender) V$sender.var else sender.var,
                         receiver.var=if(fit.vars$receiver) V$receiver.var else receiver.var,
                         sociality.var=if(fit.vars$sociality) V$sociality.var else sociality.var,
                         prior=prior,
                         opt=opt,
                         partial.eta=partial.eta,
                         up.to.const=TRUE,flyapart.penalty=flyapart.penalty)
              },
              grad.f=function(v){
                V<-unpack.optim(v,fit.vars,n,
                                p=p,
                                d=d,
                                G=G)
                # note that ergmm.lp.grad doesn't take the up.to.const parameter
                gr<-ergmm.lp.grad(Yg,response,Ym=Ym,family=family,fam.par=fam.par,
                                  beta=if(fit.vars$beta) V$beta else beta,
                                  X=X,
                                  Z=if(fit.vars$Z) V$Z else Z,
                                  sender=if(fit.vars$sender) V$sender else sender,
                                  receiver=if(fit.vars$receiver) V$receiver else receiver,
                                  sociality=if(fit.vars$sociality) V$sociality else sociality,
                                  Z.var=if(fit.vars$Z.var) V$Z.var else Z.var,
                                  Z.mean=if(fit.vars$Z.mean) V$Z.mean else Z.mean,
                                  Z.K=Z.K,
                                  sender.var=if(fit.vars$sender) V$sender.var else sender.var,
                                  receiver.var=if(fit.vars$receiver) V$receiver.var else receiver.var,
                                  sociality.var=if(fit.vars$sociality) V$sociality.var else sociality.var,
                                  prior=prior,
                                  opt=opt,
                                  partial.eta=partial.eta,
                                  flyapart.penalty=flyapart.penalty)
                pack.optim(if(fit.vars$beta) gr$beta else NULL,
                           if(fit.vars$Z) gr$Z else NULL,
                           if(fit.vars$sender) gr$sender else NULL,
                           if(fit.vars$receiver) gr$receiver else NULL,
                           if(fit.vars$sociality) gr$sociality else NULL,
                           if(fit.vars$Z.var) gr$Z.var else NULL,
                           if(fit.vars$Z.mean) gr$Z.mean else NULL,
                           if(fit.vars$sender.var) gr$sender.var else NULL,
                           if(fit.vars$receiver.var) gr$receiver.var else NULL,
                           if(fit.vars$sociality.var) gr$sociality.var else NULL)
              }
              )
           )
}

find.mle<-function(fit.vars,Yg,response=NULL,Ym=NULL,control=list(mle.maxit=200,verbose=0),
                   family="Bernoulli",fam.par=NULL,
                   beta=NULL,X=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL,
                   hessian=FALSE,mllk=TRUE,flyapart.penalty=0){
  mpe<-find.mpe(fit.vars,Yg,response=response,Ym=Ym,control=control,
                family=family,fam.par=fam.par,
                beta=beta,X=X,Z=Z,sender=sender,receiver=receiver,sociality=sociality,
                opt="lpY",
                hessian=hessian,mlp=mllk,flyapart.penalty=0)
  if(mllk)
    mpe$llk<-mpe$mlp
  mpe
}

find.mle.L<-function(model,start,given=list(),control,
                     hessian=FALSE,mllk=TRUE,Ym=NULL,flyapart.penalty=0){
  fit.vars<-list(beta=!is.null(start$beta) && is.null(given$beta),
                 Z=!is.null(start$Z) && is.null(given$Z),
                 sender=!is.null(start$sender) && is.null(given$sender),
                 receiver=!is.null(start$receiver) && is.null(given$receiver),
                 sociality=!is.null(start$sociality) && is.null(given$sociality)
                 )
  return(find.mle(fit.vars,
                  Yg=model$Yg,
                  response=model$response,
                  Ym=Ym,
                  control=control,
                  family=model$family,
                  fam.par=model$fam.par,
                  beta=(if(fit.vars$beta) start else given)$beta,
                  X=model$X,
                  Z=(if(fit.vars$Z) start else given)$Z,
                  sender=(if(fit.vars$sender) start else given)$sender,
                  receiver=(if(fit.vars$receiver) start else given)$receiver,
                  sociality=(if(fit.vars$sociality) start else given)$sociality,
                  hessian=hessian,
                  mllk=mllk,
                  flyapart.penalty=flyapart.penalty
                  )
         )
}

find.mpe<-function(fit.vars,Yg,response=NULL,Ym=NULL,control=list(mle.maxit=200,verbose=0),
                   family="Bernoulli",fam.par=NULL,
                   beta=NULL,X=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL,
                   Z.var=NULL,Z.mean=NULL,Z.K=NULL,
                   sender.var=NULL,receiver.var=NULL,sociality.var=NULL,
                   prior=NULL,opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV"),
                   hessian=FALSE,mlp=TRUE,flyapart.penalty=0){
  n<-network.size(Yg)
  p<-length(beta)
  d<-dim(Z)[2]
  G<-dim(Z.mean)[1]
  fit.vars<-reg.fit.vars(fit.vars)

  control$fnscale=-1
  control$maxit<-control$mle.maxit
  control$trace<-control$verbose
  partial.eta<-ergmm.eta(Yg,
                         beta=if(fit.vars$beta) NULL else beta,
                         X=if(fit.vars$beta) NULL else X,
                         Z=if(fit.vars$Z) NULL else Z,
                         sender=if(fit.vars$sender) NULL else sender,
                         receiver=if(fit.vars$receiver) NULL else receiver,
                         sociality=if(fit.vars$sociality) NULL else sociality) 
  
  optim.fs<-mk.lp.optim.fs(fit.vars,
                           Yg=Yg,
                           response=response,
                           Ym=Ym,
                           family=family,
                           fam.par=fam.par,
                           beta=beta,
                           X=X,
                           Z=Z,
                           sender=sender,
                           receiver=receiver,
                           sociality=sociality,
                           Z.var=Z.var,Z.mean=Z.mean,Z.K=Z.K,
                           sender.var=sender.var,receiver.var=receiver.var,sociality.var=sociality.var,
                           prior=prior,
                           opt=opt,
                           partial.eta=partial.eta,
                           flyapart.penalty=flyapart.penalty)
  
  start.vals<-pack.optim(beta,
                         Z,
                         sender,
                         receiver,
                         sociality,
                         Z.var,
                         Z.mean,
                         sender.var,
                         receiver.var,
                         sociality.var,
                         fit.vars)
  
  vmpe <- try(optim(par=start.vals,fn=optim.fs$f,gr=optim.fs$grad.f,
                    method="L-BFGS-B",
                    lower=pack.optim(
                      beta=rep(-Inf,p),
                      Z=rep(-Inf,n*d),
                      sender=rep(-Inf,n),
                      receiver=rep(-Inf,n),
                      sociality=rep(-Inf,n),
                      Z.var=rep(sqrt(.Machine$double.eps),(d>0)*max(G,1)),
                      Z.mean=rep(-Inf,d*G),
                      sender.var=sqrt(.Machine$double.eps),
                      receiver.var=sqrt(.Machine$double.eps),
                      sociality.var=sqrt(.Machine$double.eps),
                      fit.vars=fit.vars),
                    control=control,hessian=hessian))

  if(inherits(vmpe,"try-error")) return(NULL)
  mpe<-unpack.optim(vmpe$par,fit.vars,n,
                    p=p,
                    d=d,
                    G=G)

  if(!fit.vars$beta) mpe$beta<-beta
  if(!fit.vars$Z) mpe$Z<-Z
  if(!fit.vars$sender) mpe$sender<-sender
  if(!fit.vars$receiver) mpe$receiver<-receiver
  if(!fit.vars$sociality) mpe$sociality<-sociality
  if(!fit.vars$Z.var) mpe$Z.var<-Z.var
  if(!fit.vars$Z.mean) mpe$Z.mean<-Z.mean
  mpe$Z.K<-Z.K
  mpe$Z.pK<-if(!is.null(mpe$Z.K)) tabulate(Z.K)/n
  if(!fit.vars$sender.var) mpe$sender.var<-sender.var
  if(!fit.vars$receiver.var) mpe$receiver.var<-receiver.var
  if(!fit.vars$sociality.var) mpe$sociality.var<-sociality.var
  
  if(mlp)
    mpe$mlp<-with(mpe,ergmm.lp(Yg=Yg,
                               response=response,
                               Ym=Ym,
                               family=family,
                               fam.par=fam.par,
                               beta=beta,
                               X=X,
                               Z=Z,
                               sender=sender,
                               receiver=receiver,
                               sociality=sociality,
                               Z.var=Z.var,Z.mean=Z.mean,Z.K=Z.K,
                               sender.var=sender.var,receiver.var=receiver.var,sociality.var=sociality.var,
                               prior=prior,
                               opt=opt,
                               flyapart.penalty=flyapart.penalty))
  
  if(hessian) mpe$hessian<-vmpe$hessian
  class(mpe)<-"ergmm.par"
  mpe
}

find.mpe.L<-function(model,start,given=list(),prior=list(),control,fit.vars=NULL,opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV"),
                     hessian=FALSE,mlp=TRUE,Ym=NULL,flyapart.penalty=0){
  if(is.null(fit.vars)) fit.vars<-list(beta=!is.null(start$beta) && is.null(given$beta),
                                       Z=!is.null(start$Z) && is.null(given$Z),
                                       sender=!is.null(start$sender) && is.null(given$sender),
                                       receiver=!is.null(start$receiver) && is.null(given$receiver),
                                       sociality=!is.null(start$sociality) && is.null(given$sociality),
                                       Z.var=!is.null(start$Z.var) && is.null(given$Z.var),
                                       Z.mean=!is.null(start$Z.mean) && is.null(given$Z.mean),
                                       sender.var=!is.null(start$sender.var) && is.null(given$sender.var),
                                       receiver.var=!is.null(start$receiver.var) && is.null(given$receiver.var),
                                       sociality.var=!is.null(start$sociality.var) && is.null(given$sociality.var)
                                       )
  else{
    fit.vars<-reg.fit.vars(fit.vars)
    for(name in names(fit.vars)){
      if(!(name %in% names(start))) fit.vars[[name]]<-FALSE
    }
  }
              
  return(find.mpe(fit.vars,
                  Yg=model$Yg,
                  response=model$response,
                  Ym=Ym,
                  control=control,
                  family=model$family,
                  fam.par=model$fam.par,
                  beta=(if(fit.vars$beta) start else given)$beta,
                  X=model$X,
                  Z=(if(fit.vars$Z) start else given)$Z,
                  sender=(if(fit.vars$sender) start else given)$sender,
                  receiver=(if(fit.vars$receiver) start else given)$receiver,
                  sociality=(if(fit.vars$sociality) start else given)$sociality,
                  Z.var=(if(fit.vars$Z.var) start else given)$Z.var,
                  Z.mean=(if(fit.vars$Z.mean) start else given)$Z.mean,
                  Z.K=(if(!is.null(given$Z.K)) given else start)$Z.K,
                  sender.var=(if(fit.vars$sender.var) start else given)$sender.var,
                  receiver.var=(if(fit.vars$receiver.var) start else given)$receiver.var,
                  sociality.var=(if(fit.vars$sociality.var) start else given)$sociality.var,
                  prior=prior,
                  opt=opt,
                  hessian=hessian,
                  mlp=mlp,
                  flyapart.penalty=flyapart.penalty
                  )
         )
}

ergmm.loglike.L<-function(model,theta,Ym=NULL,flyapart.penalty=0,up.to.const=FALSE)
  ergmm.loglike(Yg=model$Yg,response=model$response,Ym=Ym,
                family=model$family,fam.par=model$fam.par,
                beta=theta$beta,X=model$X,Z=theta$Z,
                sender=theta$sender,receiver=theta$receiver,sociality=theta$sociality,
                flyapart.penalty=flyapart.penalty,up.to.const=up.to.const)

ergmm.loglike.grad.L<-function(model,theta,Ym=NULL,flyapart.penalty=0)
  ergmm.loglike.grad(Yg=model$Yg,response=model$response,Ym=Ym,
                     family=model$family,fam.par=model$fam.par,
                     beta=theta$beta,X=model$X,Z=theta$Z,
                     sender=theta$sender,receiver=theta$receiver,sociality=theta$sociality,
                     flyapart.penalty=flyapart.penalty)


ergmm.lp<-function(Yg,response=NULL,Ym=NULL,family="Bernoulli",fam.par=NULL,beta=NULL,X=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL,
                   Z.var=NULL,Z.mean=NULL,Z.K=NULL,
                   sender.var=NULL,receiver.var=NULL,sociality.var=NULL,
                   prior,opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV"),
                   partial.eta=NULL,up.to.const=FALSE,flyapart.penalty=0){
  lpY<-if("lpY" %in% opt) ergmm.loglike(Yg,response=response,Ym=Ym,family=family,fam.par=fam.par,
                                        beta=beta,X=X,Z=Z,sender=sender,receiver=receiver,sociality=sociality,
                                        partial.eta=partial.eta,up.to.const=up.to.const,flyapart.penalty=0) else 0
  
  lpZ<-if("lpZ" %in% opt) ergmm.lpZ(Z=Z,Z.var=Z.var,Z.mean=Z.mean,Z.K=Z.K) else 0
  lpRE<-if("lpRE" %in% opt) ergmm.lpRE(sender=sender,receiver=receiver,sociality=sociality,
                                       sender.var=sender.var,receiver.var=receiver.var,sociality.var=sociality.var) else 0
  lpBeta<-if("lpBeta" %in% opt) ergmm.lpBeta(beta,prior$beta.mean,prior$beta.var) else 0
  lpREV<-if("lpREV" %in% opt) ergmm.lpREV(sender.var=sender.var,receiver.var=receiver.var,sociality.var=sociality.var,
                                          sender.var.prior=prior$sender.var,receiver.var.prior=prior$receiver.var,sociality.var.prior=prior$sociality.var,
                                          sender.var.df=prior$sender.var.df,receiver.var.df=prior$receiver.var.df,sociality.var.df=prior$sociality.var.df) else 0
  
  lpLV<-if("lpLV" %in% opt) ergmm.lpLV(Z.var=Z.var,Z.mean=Z.mean,Z.var.prior=prior$Z.var,Z.var.df=prior$Z.var.df,Z.mean.var=prior$Z.mean.var) else 0
  
  lpAll<-lpY+lpZ+lpRE+lpBeta+lpREV+lpLV

  if(flyapart.penalty) return(lpAll+flyapart.penalty*(sum(beta^2)+
                                                      sum(Z^2)+
                                                      sum(sender^2)+
                                                      sum(receiver^2)+
                                                      sum(sociality^2)+
                                                      sum(Z.var^2)+
                                                      sum(Z.mean^2)+
                                                      sum(sender.var^2)+
                                                      sum(receiver.var^2)+
                                                      sum(sociality.var^2)))
                                                      
  return(lpAll)
}

ergmm.lp.L<-function(model,theta,prior,Ym=NULL,opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV"),up.to.const=FALSE,flyapart.penalty=0){
  ergmm.lp(Yg=model$Yg,response=model$response,Ym=Ym,
               family=model$family,fam.par=model$fam.par,
               beta=theta$beta,X=model$X,Z=theta$Z,
               sender=theta$sender,receiver=theta$receiver,sociality=theta$sociality,
               Z.var=theta$Z.var,Z.mean=theta$Z.mean,Z.K=theta$Z.K,
               sender.var=theta$sender.var,receiver.var=theta$receiver.var,sociality.var=theta$sociality.var,
               prior,opt=opt,
               up.to.const=up.to.const,flyapart.penalty=flyapart.penalty)
}

add.lists<-function(...){
  out<-list()
  for(l in list(...)){
    for(name in names(l)){
      if(name %in% names(out))
        out[[name]]<-out[[name]]+l[[name]]
      else
        out[[name]]<-l[[name]]
    }
  }
  out
}

ergmm.lp.grad<-function(Yg,response=NULL,Ym=NULL,family="Bernoulli",fam.par=NULL,beta=NULL,X=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL,
                        Z.var=NULL,Z.mean=NULL,Z.K=NULL,
                        sender.var=NULL,receiver.var=NULL,sociality.var=NULL,
                        prior,opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV"),
                        partial.eta=NULL,flyapart.penalty=0){
  
  
  grad<-add.lists(if("lpY" %in% opt) if(!is.null(beta)||!is.null(Z)||!is.null(sender)||!is.null(receiver)||!is.null(sociality)) ergmm.loglike.grad(Yg,response=response,Ym=Ym,family=family,fam.par=fam.par,
                                                                                                                                beta=beta,X=X,Z=Z,sender=sender,receiver=receiver,sociality=sociality,
                                                                                                                                partial.eta=partial.eta,flyapart.penalty=0),
                  if("lpZ" %in% opt) ergmm.lpZ.grad(Z=Z,Z.var=Z.var,Z.mean=Z.mean,Z.K=Z.K),
                  if("lpRE" %in% opt) ergmm.lpRE.grad(sender=sender,receiver=receiver,sociality=sociality,
                                  sender.var=sender.var,receiver.var=receiver.var,sociality.var=sociality.var),
                  if("lpBeta" %in% opt) ergmm.lpBeta.grad(beta,prior$beta.mean,prior$beta.var),
                  if("lpREV" %in% opt) ergmm.lpREV.grad(sender.var=sender.var,receiver.var=receiver.var,sociality.var=sociality.var,
                                   sender.var.prior=prior$sender.var,receiver.var.prior=prior$receiver.var,sociality.var.prior=prior$sociality.var,
                                   sender.var.df=prior$sender.var.df,receiver.var.df=prior$receiver.var.df,sociality.var.df=prior$sociality.var.df),
                  if("lpLV" %in% opt) ergmm.lpLV.grad(Z.var=Z.var,Z.mean=Z.mean,Z.var.prior=prior$Z.var,Z.var.df=prior$Z.var.df,Z.mean.var=prior$Z.mean.var))
                  
  if(!is.null(grad$beta))grad$beta<-grad$beta-flyapart.penalty*2*beta
  if(!is.null(grad$Z))grad$Z<-grad$Z-flyapart.penalty*2*Z
  if(!is.null(grad$sender))grad$sender<-grad$sender-flyapart.penalty*2*sender
  if(!is.null(grad$receiver))grad$receiver<-grad$receiver-flyapart.penalty*2*receiver
  if(!is.null(grad$sociality))grad$sociality<-grad$sociality-flyapart.penalty*2*sociality
  if(!is.null(grad$Z.var))grad$Z.var<-grad$Z.var-flyapart.penalty*2*Z.var
  if(!is.null(grad$Z.mean)) grad$Z.mean<-grad$Z.mean-flyapart.penalty*2*Z.mean
  if(!is.null(grad$sender.var))grad$sender.var<-grad$sender.var-flyapart.penalty*2*sender.var
  if(!is.null(grad$receiver.var))grad$receiver.var<-grad$receiver.var-flyapart.penalty*2*receiver.var
  if(!is.null(grad$sociality.var))grad$sociality.var<-grad$sociality.var-flyapart.penalty*2*sociality.var
  grad
}

ergmm.lp.grad.L<-function(model,theta,prior,Ym=NULL,opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV"),flyapart.penalty=0){
  ergmm.lp.grad(Yg=model$Yg,response=model$response,Ym=Ym,
                    family=model$family,fam.par=model$fam.par,
                    beta=theta$beta,X=model$X,Z=theta$Z,
                    sender=theta$sender,receiver=theta$receiver,sociality=theta$sociality,
                    Z.var=theta$Z.var,Z.mean=theta$Z.mean,Z.K=theta$Z.K,
                    sender.var=theta$sender.var,receiver.var=theta$receiver.var,sociality.var=theta$sociality.var,
                    prior,opt=opt,
                    flyapart.penalty=flyapart.penalty)
}

ergmm.lp.grad.approx<-function(Yg,response=NULL,Ym=NULL,family="Bernoulli",fam.par=NULL,beta=NULL,X=NULL,Z=NULL,sender=NULL,receiver=NULL,sociality=NULL,
                               Z.var=NULL,Z.mean=NULL,Z.K=NULL,
                               sender.var=NULL,receiver.var=NULL,sociality.var=NULL,
                               prior,opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV"),
                               partial.eta=NULL,flyapart.penalty=0,delta=sqrt(.Machine$double.eps),which.vars=FIT_ALL){
  which.vars<-reg.fit.vars(which.vars)
  fit.vars<-list(beta=!is.null(beta),
                 Z=!is.null(Z),
                 sender=!is.null(sender),
                 receiver=!is.null(receiver),
                 sociality=!is.null(sociality),
                 Z.var=!is.null(Z.var),
                 Z.mean=!is.null(Z.mean),
                 sender.var=!is.null(sender.var),
                 receiver.var=!is.null(receiver.var),
                 sociality.var=!is.null(sociality.var))
  which.vars.v<-pack.optim(beta=rep(which.vars$beta,length(beta)),
                           Z=rep(which.vars$Z,length(Z)),
                           sender=rep(which.vars$sender,length(sender)),
                           receiver=rep(which.vars$receiver,length(receiver)),
                           sociality=rep(which.vars$sociality,length(sociality)),
                           Z.var=rep(which.vars$Z.var,length(Z.var)),
                           Z.mean=rep(which.vars$Z.mean,length(Z.mean)),
                         sender.var=rep(which.vars$sender.var,length(sender.var)),
                           receiver.var=rep(which.vars$receiver.var,length(receiver.var)),
                           sociality.var=rep(which.vars$sociality.var,length(sociality.var)),
                           fit.vars)
  n<-network.size(Yg)
  p<-length(beta)
  d<- if(is.null(Z)) 0 else dim(Z)[2]
  G<-if(is.null(Z.mean)) 0 else dim(Z.mean)[1]

  model<-list(Yg=Yg,response=response,Ym=Ym,family=family,fam.par=fam.par,X=X,
              sender=!is.null(sender)||!is.null(sender.var),
              receiver=!is.null(receiver)||!is.null(receiver.var),
              sociality=!is.null(sociality)||!is.null(receiver.var))
  
  fit.vars<-reg.fit.vars(fit.vars)
  v<-pack.optim(beta,Z,sender,receiver,sociality,Z.var,Z.mean,sender.var,receiver.var,sociality.var,fit.vars)
  dlpdv<-numeric(length(v))
  
  for(i in 1:length(v)){
    if(!which.vars.v[i]) next
    v.m<-v.p<-v
    v.m[i]<-v.m[i]-delta
    v.p[i]<-v.p[i]+delta
    lp.m<-ergmm.lp.L(model,c(unpack.optim(v.m,fit.vars,n,p,d,G),list(Z.K=Z.K)),prior,Ym=Ym,opt=opt,flyapart.penalty=flyapart.penalty)
    lp.p<-ergmm.lp.L(model,c(unpack.optim(v.p,fit.vars,n,p,d,G),list(Z.K=Z.K)),prior,Ym=Ym,opt=opt,flyapart.penalty=flyapart.penalty)
    dlpdv[i]<-(lp.p-lp.m)/(2*delta)
  }
  
  return(unpack.optim(dlpdv,fit.vars,n,p,d,G))
}


ergmm.lp.grad.approx.L<-function(model,theta,prior,Ym=NULL,opt=c("lpY","lpZ","lpBeta","lpRE","lpREV","lpLV"),flyapart.penalty=0,delta=sqrt(.Machine$double.eps),which.vars=FIT_ALL){
  ergmm.lp.grad.approx(Yg=model$Yg,response=model$response,Ym=Ym,
                       family=model$family,fam.par=model$fam.par,
                       beta=theta$beta,X=model$X,Z=theta$Z,
                       sender=theta$sender,receiver=theta$receiver,sociality=theta$sociality,
                       Z.var=theta$Z.var,Z.mean=theta$Z.mean,Z.K=theta$Z.K,
                       sender.var=theta$sender.var,receiver.var=theta$receiver.var,sociality.var=theta$sociality.var,
                       prior,opt=opt,
                       flyapart.penalty=flyapart.penalty,delta=delta,which.vars=which.vars)
}

ergmm.lpZ<-function(Z=NULL,Z.var=NULL,Z.mean=NULL,Z.K=NULL){
  if(is.null(Z)||is.null(Z.var)) return(0)
  n<-dim(Z)[1]
  d<-dim(Z)[2]
  if(is.null(Z.K)){
    if(!is.null(Z.mean)) stop("Given cluster means without cluster assignments!")
    Z.K<-rep(1,n)
    Z.mean<-matrix(0,nrow=1,ncol=d)
  }
  sum(dnorm(Z,Z.mean[Z.K,],matrix(sqrt(Z.var[Z.K]),nrow=n,ncol=d,byrow=FALSE),TRUE))
}

ergmm.lpZ.grad<-function(Z=NULL,Z.var=NULL,Z.mean=NULL,Z.K=NULL){
  if(is.null(Z)||is.null(Z.var)) return(list())
  n<-dim(Z)[1]
  d<-dim(Z)[2]
  G<-if(is.null(Z.K)) 0 else dim(Z.mean)[1]
  deriv<-list()
  if(is.null(Z.K)){
    if(!is.null(Z.mean)) stop("Given cluster means without cluster assignments!")
    Z.K<-rep(1,n)
    Z.mean<-matrix(0,nrow=1,ncol=d)
  }
  Z.dev<-(Z-Z.mean[Z.K,])/matrix(Z.var[Z.K],nrow=n,ncol=d,byrow=FALSE)
  deriv$Z<--Z.dev
  deriv$Z.var<-sapply(1:max(G,1),
                        function(g)
                        (sum(Z.dev[Z.K==g,,drop=FALSE]^2)-d*sum(Z.K==g)/Z.var[g])/2)
  
  if(G) deriv$Z.mean<-t(sapply(1:G,function(g) apply(Z.dev[Z.K==g,,drop=FALSE],2,sum)))

  deriv
}

ergmm.lpZ.L<-function(theta){
  ergmm.lpZ(theta$Z,theta$Z.var,theta$Z.mean,theta$Z.K)
}

ergmm.lpZ.grad.L<-function(theta){
  ergmm.lpZ.grad(theta$Z,theta$Z.var,theta$Z.mean,theta$Z.K)
}

ergmm.lpRE<-function(sender=NULL,receiver=NULL,sociality=NULL,sender.var=NULL,
                     receiver.var=NULL,sociality.var=NULL){
  ({if(is.null(sender) || is.null(sender.var)) 0 else sum(dnorm(sender,0,sqrt(sender.var),TRUE))}+
   {if(is.null(receiver) || is.null(receiver.var)) 0 else sum(dnorm(receiver,0,sqrt(receiver.var),TRUE))}+
   {if(is.null(sociality) || is.null(sociality.var)) 0 else sum(dnorm(sociality,0,sqrt(sociality.var),TRUE))})
}

ergmm.lpRE.grad<-function(sender=NULL,receiver=NULL,sociality=NULL,sender.var=NULL,
                          receiver.var=NULL,sociality.var=NULL){
  deriv<-list()
  if(!is.null(sender) && !is.null(sender.var)){
    deriv$sender<--sender/sender.var
    deriv$sender.var<-(sum(sender^2)/sender.var-length(sender))/sender.var/2
  }
  if(!is.null(receiver) && !is.null(receiver.var)){
    deriv$receiver<--receiver/receiver.var
    deriv$receiver.var<-(sum(receiver^2)/receiver.var-length(receiver))/receiver.var/2
  }
  if(!is.null(sociality) && !is.null(sociality.var)){
    deriv$sociality<--sociality/sociality.var
    deriv$sociality.var<-(sum(sociality^2)/sociality.var-length(sociality))/sociality.var/2
  }
  
  deriv
}

ergmm.lpRE.L<-function(theta){
  ergmm.lpRE(theta$sender,theta$receiver,theta$sociality,
             theta$sender.var,theta$receiver.var,theta$sociality.var)
}

ergmm.lpRE.grad.L<-function(theta){
  ergmm.lpRE.grad(theta$sender,theta$receiver,theta$sociality,
                  theta$sender.var,theta$receiver.var,theta$sociality.var)
}

ergmm.lpBeta<-function(beta=NULL,beta.mean=NULL,beta.var=NULL){
  if(is.null(beta)||is.null(beta.mean)||is.null(beta.var)) return(0)
  sum(dnorm(beta,beta.mean,sqrt(beta.var),TRUE))
}

ergmm.lpBeta.grad<-function(beta=NULL,beta.mean=NULL,beta.var=NULL){
  if(is.null(beta)||is.null(beta.mean)||is.null(beta.var)) return(list())
  deriv<-list(beta=-(beta-beta.mean)/beta.var)
  deriv
}
ergmm.lpBeta.L<-function(theta,prior){
  ergmm.lpBeta(theta$beta,prior$beta.mean,prior$beta.var)
}
ergmm.lpBeta.grad.L<-function(theta,prior){
  ergmm.lpBeta.grad(theta$beta,prior$beta.mean,prior$beta.var)
}


dsclinvchisq<-function(x,df,scale=1,log=FALSE){
  if(log) dchisq(df*scale/x,df,log=TRUE)+log(df)+log(scale)-2*log(x)
  else dchisq(df*scale/x,df,log=FALSE)*df*scale/x/x
}

ergmm.lpREV<-function(sender.var=NULL,receiver.var=NULL,sociality.var=NULL,
                      sender.var.prior=NULL,receiver.var.prior=NULL,sociality.var.prior=NULL,
                      sender.var.df=NULL,receiver.var.df=NULL,sociality.var.df=NULL){
  ({if(is.null(sender.var)) 0 else dsclinvchisq(sender.var,sender.var.df,sender.var.prior,TRUE)}+
   {if(is.null(receiver.var)) 0 else dsclinvchisq(receiver.var,receiver.var.df,receiver.var.prior,TRUE)}+
   {if(is.null(sociality.var)) 0 else dsclinvchisq(sociality.var,sociality.var.df,sociality.var.prior,TRUE)})
}

ergmm.lpREV.grad<-function(sender.var=NULL,receiver.var=NULL,sociality.var=NULL,
                           sender.var.prior=NULL,receiver.var.prior=NULL,sociality.var.prior=NULL,
                           sender.var.df=NULL,receiver.var.df=NULL,sociality.var.df=NULL){
  deriv<-list()
  if(!is.null(sender.var)) deriv$sender.var<-sender.var.df*sender.var.prior/sender.var^2/2-(sender.var.df/2+1)/sender.var
  if(!is.null(receiver.var)) deriv$receiver.var<-receiver.var.df*receiver.var.prior/receiver.var^2/2-(receiver.var.df/2+1)/receiver.var
  if(!is.null(sociality.var)) deriv$sociality.var<-sociality.var.df*sociality.var.prior/sociality.var^2/2-(sociality.var.df/2+1)/sociality.var
  
  deriv
}

ergmm.lpREV.L<-function(theta,prior){
  ergmm.lpREV(theta$sender.var,theta$receiver.var,theta$sociality.var,
              prior$sender.var,prior$receiver.var,prior$sociality.var,
              prior$sender.var.df,prior$receiver.var.df,prior$sociality.var.df)
}
ergmm.lpREV.grad.L<-function(theta,prior){
  ergmm.lpREV.grad(theta$sender.var,theta$receiver.var,theta$sociality.var,
                   prior$sender.var,prior$receiver.var,prior$sociality.var,
                   prior$sender.var.df,prior$receiver.var.df,prior$sociality.var.df)
}

ergmm.lpLV<-function(Z.var=NULL,Z.mean=NULL,Z.var.prior=NULL,Z.var.df=NULL,Z.mean.var=NULL){
  (if(is.null(Z.var)) 0 else sum(dsclinvchisq(Z.var,Z.var.df,Z.var.prior,log=TRUE))+
   {if(is.null(Z.mean)) 0 else sum(dnorm(Z.mean,0,sqrt(Z.mean.var),log=TRUE))})
}
ergmm.lpLV.L<-function(theta,prior){
  ergmm.lpLV(theta$Z.var,theta$Z.mean,prior$Z.var,prior$Z.var.df,prior$Z.mean.var)
}
ergmm.lpLV.grad<-function(Z.var=NULL,Z.mean=NULL,Z.var.prior=NULL,Z.var.df=NULL,Z.mean.var=NULL){
  deriv<-list()
  if(!is.null(Z.var)) deriv$Z.var<-Z.var.df*Z.var.prior/Z.var^2/2-(Z.var.df/2+1)/Z.var
  if(!is.null(Z.mean)) deriv$Z.mean<--Z.mean/Z.mean.var
  deriv
}
ergmm.lpLV.grad.L<-function(theta,prior){
  ergmm.lpLV.grad(theta$Z.var,theta$Z.mean,prior$Z.var,prior$Z.var.df,prior$Z.mean.var)
}
