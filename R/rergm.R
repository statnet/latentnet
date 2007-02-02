if(!exists("rergm", mode="function")){
  rergm <- function(object, ...)
    UseMethod("rergm")
}

rergm.ergmm<-function(ergmm.fit,n=1){
  l<-list()
  for(i in 1:n){
    iter<-floor(runif(1,1,ergmm.fit$control$samplesize+1))
    l[[i]]<-rergm.1(ergmm.fit$model,ergmm.fit$samples[[iter]],ergmm.fit$prior)
  }
  if(n==1) return(l[[1]])
  else{
    attr(l,"class")<-"network.series"
    return(l)
  }
}

rergm.list<-function(model,par,prior=ergmm.prior(),n=1){
  l<-list()
  for(i in 1:n){
    l[[i]]<-rergm.1(model,par,prior)
  }
  if(n==1) return(l[[1]])
  else{
    attr(l,"class")<-"network.series"
    return(l)
  }
}
  
rergm.1<-function(model,par,prior=ergmm.prior()){
  nv<-model$Yg$gal$n
  mypar<-par
  
  if(length(model$X)>0 && is.null(mypar$beta))
    mypar$beta<-rnorm(length(model$X),prior$beta.mean,sqrt(prior$beta.var))

  if(model$d>0 && is.null(mypar$Z)){
    if(model$G>0){
      if(is.null(mypar$Z.mean))
        mypar$Z.mean<-matrix(rnorm(model$G*model$d,0,sqrt(prior$Z.mean.var)),nrow=model$G)
      if(is.null(mypar$Z.K))
        mypar$Z.K<-sample(1:model$G,nv,replace=TRUE)
    }
    
    if(is.null(mypar$Z.var))
      mypar$Z.var<-prior$Z.var*prior$Z.var.df/rchisq(max(model$G,1),prior$Z.var.df)

    mypar$Z<-matrix(rnorm(nv*model$d,
                            if(model$G>0) mypar$Z.mean[mypar$Z.K,] else 0,
                            if(model$G>0) mypar$Z.var[mypar$Z.K,] else mypar$Z.var
                            ),nrow=nv)
  }

  if(model$sociality && is.null(mypar$sociality)){
    if(is.null(mypar$sociality.var))
      mypar$sociality.var<-with(prior,sociality.var*sociality.var.df/rchisq(1,sociality.var.df))
    model$sociality<-rnorm(nv,0,sqrt(mypar$sociality.var))
  }

  if(model$sender && is.null(mypar$sender)){
    if(is.null(mypar$sender.var))
      mypar$sender.var<-with(prior,sender.var*sender.var.df/rchisq(1,sender.var.df))
    mypar$sender<-rnorm(nv,0,sqrt(mypar$sender.var))
  }
  
  if(model$receiver && is.null(mypar$receiver)){
    if(is.null(mypar$receiver.var))
      mypar$receiver.var<-with(prior,receiver.var*receiver.var.df/rchisq(1,receiver.var.df))
    mypar$receiver<-rnorm(nv,0,sqrt(mypar$receiver.var))
  }

  eta<-ergmm.eta(model$Yg,
                 beta=mypar$beta,
                 X=model$X,
                 Z=mypar$Z,
                 sender=mypar$sender,
                 receiver=mypar$receiver,
                 sociality=mypar$sociality)

  sm<-mk.rsm.f(model$family)(eta,model$fam.par)
  if(!has.loops(model$Yg))
    diag(sm)<-0
  net<-with(model,as.network.matrix(sm,matrix.type="adjacency",
                                    directed=is.directed(Yg),
                                    loops=has.loops(Yg)))
  if(!is.null(model$response))
    net<-set.edge.value(net,model$response,sm)

  attr(net,"ergmm.par")<-mypar

  net
}
