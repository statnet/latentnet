# Utilities for dealing with MCMC output produced by *.MCMC.C functions.

ERGMM.PAR_VAR_NAMES<-c("beta","Z","sender","receiver","sociality",
                       "Z.var","Z.mean","Z.K",
                       "sender.var","receiver.var","sociality.var")
ERGMM.PAR_LLK_NAMES<-c("beta","Z","sender","receiver","sociality")

del.iteration<-function(mcmcsamples,i){
  for(name in names(mcmcsamples)){
    if(length(mcmcsamples[[name]])>0){
      if(length(dim(mcmcsamples[[name]]))<=1) mcmcsamples[[name]]<-mcmcsamples[[name]][-i]
      else if(length(dim(mcmcsamples[[name]]))==2) mcmcsamples[[name]]<-mcmcsamples[[name]][-i,,drop=FALSE]
      else if(length(dim(mcmcsamples[[name]]))==3) mcmcsamples[[name]]<-mcmcsamples[[name]][-i,,,drop=FALSE]
    }
  }
  mcmcsamples
}

seldrop<-function(x,i){
  array(c(x),dim=dim(x)[-i])
}

ergmm.par.blank<-function(){
  x<-list()
  class(x)<-"ergmm.par"
  x
}

as.ergmm.par.list<-function(x){
  class(x)<-"ergmm.par"
  x
}

"[.ergmm.par.list"<-function(x,i){
  if(!is.numeric(i)) stop("Index vector to the '[' operator of an ergmm.par.list must be of mode numeric or integer.")
  l<-list()
  
  for(name in names(x)){
    if(length(x[[name]])>0){
      d<-dim(x[[name]])
      if(length(d)<=1) l[[name]]<-x[[name]][i]
      else if(length(d)==2) l[[name]]<-x[[name]][i,,drop=FALSE]
      else if(length(d)==3) l[[name]]<-x[[name]][i,,,drop=FALSE]
    }
  }
  class(l)<-"ergmm.par.list"
  l
}

length.ergmm.par.list<-function(x){
  if(is.null(dim(x[[names(x)[1]]]))) length(x[[names(x)[1]]])
  else dim(x[[names(x)[1]]])[1]
}

"[[.ergmm.par.list"<-"$.ergmm.par.list"<-function(x,i){
  ## Delete its class, to keep it from recursing.
  tmp<-class(x)
  class(x)<-NULL
  if(class(i)=="character"){ ## If the index is a character, return all the samples for the corresponding variable.
    xi<-x[i][[1]]
    class(x)<-tmp
    return(xi)
  }
  else{  ## If it's a number, return a configuration with that iteration number.
    l<-list()
    ## Do NOT seldrop 1D parameters. Those actually need to become vectors.
    ## (As opposed to becoming 1*n matrices.)
    
    for(name in names(x)){
      if(length(x[[name]])>0){
        if(length(dim(x[[name]]))<=1) l[[name]]<-x[[name]][i]
        else if(length(dim(x[[name]]))==2) l[[name]]<-x[[name]][i,]
        else if(length(dim(x[[name]]))==3) l[[name]]<-seldrop(x[[name]][i,,,drop=FALSE],1)
      }
    }
    class(l)<-"ergmm.par"
    class(x)<-tmp
    return(l)
  }
}

"[[.ergmm.par"<-"$.ergmm.par"<-function(x,i){
  x[i][[1]]
}

stack.ergmm.par.list.list<-function(mcmcList){
  require(abind)
  mcmcsamples<-list()

  for(name in names(mcmcList[[1]]))
    mcmcsamples[[name]]<-abind(sapply(1:length(mcmcList),
                                      function(i) mcmcList[[i]][[name]],
                                      simplify=FALSE),along=1)

  attr(mcmcsamples,"breaks")<-cumsum(c(sapply(1:(length(mcmcList)
                                                 ),
                                              function(i) length(mcmcList[[i]]),
                                              simplify=FALSE)))
  class(mcmcsamples)<-"ergmm.par.list"
  mcmcsamples
}

unstack.ergmm.par.list<-function(mcmcsamples){
  mcmcList<-list()

  if(is.null(attr(mcmcsamples,"breaks"))){
    mcmcList[[1]]<-mcmcsamples
  }
  else{  
    breaks<-c(0,attr(mcmcsamples,"breaks"))
    
    for(i in 1:length(breaks[-1])){
      mcmcList[[i]]<-mcmcsamples[(breaks[i]+1):breaks[i+1]]
    }
  }
  mcmcList
}

as.mcmc.list.ergmm.par.list<-function(x,which.vars,start=1,thin=1){
  require(coda)
  x<-unstack(x)
  m.l<-list()
  for(thread in 1:length(x)){
    S<-length(x[[thread]])
    m<-matrix(numeric(0),S,0)
    for(name in names(which.vars)){
      if(is.null(x[[thread]][[name]])) next
      if(length(dim(x[[thread]][[name]]))<=1){
        m2<-cbind(x[[thread]][[name]])
        colnames(m2)<-name
        m<-cbind(m,m2)
      }else if(length(dim(x[[thread]][[name]]))==2){
        m2<-x[[thread]][[name]][,c(which.vars[[name]]),drop=FALSE]
        colnames(m2)<-paste(name,sapply(which.vars[[name]],function(x) paste(x,sep='.')),sep='.')
        m<-cbind(m,m2)
      }else if(length(dim(x[[thread]][[name]]))==3){
        for(i in 1:dim(which.vars[[name]])[1]){
          i<-which.vars[[name]][i,]
          m2<-cbind(x[[thread]][[name]][,i[1],i[2]])
          colnames(m2)<-paste(name,paste(i,sep='.',collapse='.'),sep='.')
          m<-cbind(m,m2)
        }
      }
    }
    m<-mcmc(m,start=start,thin=thin)
    m.l[[thread]]<-m
  }
  eval(as.call(c(mcmc.list,m.l)))
}
