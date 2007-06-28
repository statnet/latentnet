ergmm.build.model <- function(formula,response,family,fam.par,orthogonalize,prior){
  
  terms<-terms(formula)

  if(!attr(terms,"response") || terms[[1]]!="~") stop("Formula must be of form 'network ~ model'.")

  Yg <- try(as.network(eval(terms[[2]],attr(terms,".Environment"))))
  if(inherits(Yg,"try-error")){
    stop("Invalid network. Is the left-hand-side of the formula correct?")
  }

  model<-list(formula=formula,
              Yg=Yg,
              response=response,
              family=family,
              fam.par=fam.par,
              coef.names=character(0),
              X=list(),
              p=0,
              d=0,
              G=0,
              sender=FALSE,
              receiver=FALSE,
              sociality=FALSE,
              intercept=as.logical(attr(terms,"intercept")),
              prior=list() ## Only here for convenience.
              )



  if(model$intercept){
    model<-InitErgmm.latentcov(model,observed.dyads(Yg),"density")
  }
              
  for (term in as.list(attr(terms,"variables"))[-(1:2)]){
    if (is.call(term)){
      init.call<-list()
      init.call<-list(as.name(paste("InitErgmm.", term[[1]], sep = "")),
                      model=model)
      
      init.call<-c(init.call,as.list(term)[-1])
    }else{
      init.call <- list(as.name(paste("InitErgmm.", term, sep = "")),model=model)
    }
    model <- eval(as.call(init.call), attr(terms,".Environment"))
  }
  
  if(orthogonalize && model$p>1)
    model$X<-GS.orth.matrix(model$X,observed.dyads(Yg))

  for(name in names(prior)){
    model$prior[[name]]<-prior[[name]]
  }

  prior<-model$prior
  model$prior<-NULL

  if(!("Z.var" %in% names(prior))) prior$Z.var<-prior$Z.var.mul*(network.size(model$Yg)/max(1,model$G))^(2/model$d)
  if(!("Z.mean.var" %in% names(prior))) prior$Z.mean.var<-prior$Z.mean.var.mul*prior$Z.var*max(1,model$G)^(2/model$d)
  if(!("Z.var.df" %in% names(prior))) prior$Z.var.df<-prior$Z.var.df.mul*sqrt(network.size(model$Yg)/max(1,model$G))
  if(!("Z.pK" %in% names(prior))) prior$Z.pK<-prior$Z.pK.mul*sqrt(network.size(model$Yg)/max(1,model$G))

  class(model)<-"ergmm.model"  
  list(model=model,prior=prior)
}

GS.orth.matrix<-function(X,missing){
  p<-length(X)
  if(p<=1) return(X)
  
  Y<-list()
  Y[[1]]<-X[[1]]
  for(i in 2:p){
    Y[[i]]<-X[[i]]
    
    for(j in 1:(i-1)){
      Yv<-Y[[j]][!missing]
      Xv<-Y[[i]][!missing]
      Y[[i]]<-Y[[i]]-crossprod(Yv,Xv)/sqrt(crossprod(Xv))*X[[j]]
    }
  }
  Y
}
