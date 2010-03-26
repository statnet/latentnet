ergmm.get.model <- function(formula,response,family,fam.par,prior){
  
  terms<-terms(formula)

  if(!attr(terms,"response") || terms[[1]]!="~") stop("Formula must be of form 'network ~ model'.")

  Yg <- try(as.network(eval(terms[[2]],attr(terms,".Environment"))))
  if(inherits(Yg,"try-error")){
    stop("Invalid network. Is the left-hand-side of the formula correct?")
  }

  model<-list(formula=formula,
              Yg=Yg,
              Ym=getYm(Yg,response),
              response=response,
              family=family,
              familyID=family.IDs[[family]],
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

  model<-fam.par.check(model)
  
  if(model[["intercept"]]){
    model<-InitErgmm.latentcov(model,matrix(1,network.size(Yg),network.size(Yg)),"edges")
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

  if(model[["d"]]>0) model[["latentID"]] <- latent.effect.IDs[[model[["latent"]]]]

  if(prior[["adjust.beta.var"]]) model[["prior"]][["beta.var"]]<-model[["prior"]][["beta.var"]]/sapply(1:model[["p"]],function(i) mean((model[["X"]][[i]][observed.dyads(model[["Yg"]])])^2))
  
  for(name in names(prior)){
    model[["prior"]][[name]]<-prior[[name]]
  }

  prior<-model[["prior"]]
  model[["prior"]]<-NULL

  beta.eff<-get.beta.eff(model)
  for(re in names(beta.eff))
    model[[paste("beta.eff",re,sep=".")]]<-beta.eff[[re]]
  
  class(model)<-"ergmm.model"  
  list(model=model,prior=prior)
}

get.beta.eff<-function(model){
  n<-network.size(model[["Yg"]])
  out<-list(sender = if(model[["sender"]]) t(sapply(1:model[["p"]],function(k) apply(model[["X"]][[k]],1,mean))),
            receiver = if(model[["receiver"]]) t(sapply(1:model[["p"]],function(k) apply(model[["X"]][[k]],2,mean))),
            sociality = if(model[["sociality"]]) t(sapply(1:model[["p"]],function(k) apply(model[["X"]][[k]],1,mean))))
  for(re in names(out))
    if(is.null(out[[re]])) out[[re]]<-NULL
    else if(model[["p"]]>1)
      for(k1 in 2:model[["p"]])
        for(k2 in 1:(k1-1)){
          utu<-crossprod(out[[re]][k2,],out[[re]][k2,])
          if(isTRUE(all.equal(utu,0))) break;
          out[[re]][k1,]<-out[[re]][k1,]-crossprod(out[[re]][k2,],out[[re]][k1,])/utu*out[[re]][k2,]
        }
  ## Rows of 0s will break the tuner, and don't actually help, so, we drop them.
  for(re in names(out)){
    no.eff<-apply(out[[re]],1,function(x) isTRUE(all.equal(x,rep(0,n))))
    out[[re]]<-out[[re]][!no.eff,,drop=FALSE]
  }
  
  out
}
