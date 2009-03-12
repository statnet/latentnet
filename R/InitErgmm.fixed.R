InitErgmm.latentcov<-function (model, x, attrname=NULL,
                               mean=0, var=9) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("latentcov() model term expected between 1 and 4 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  #Coerce x to an adjacency matrix
  if(is.network(x)){
    xm<-as.matrix.network(x,matrix.type="adjacency",attrname)
    cn<-attrname
  }else if(is.character(x)){
    xm<-as.matrix.network(model[["Yg"]],matrix.type="adjacency",x)
    cn<-x
  }else{
    xm<-as.matrix(x)
    cn<-if(!is.null(attrname)) attrname else paste("matrix",length(model[["X"]])+1)
  }

  model[["p"]]<-model[["p"]]+1
  model[["X"]]<-c(model[["X"]],list(xm))
  model[["coef.names"]]<-c(model[["coef.names"]], cn)

  model[["prior"]][["beta.mean"]]<-c(model[["prior"]][["beta.mean"]],mean)
  model[["prior"]][["beta.var"]]<-c(model[["prior"]][["beta.var"]],var)

  model
}
InitErgmm.absdiff<-function (model, attrname,
                             mean=0, var=9) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:4))
    stop(paste("absdiff() model term expected between 1 and 3 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  x<-model[["Yg"]] %v% attrname
  if(!is.vector(x)||!is.numeric(x))
    stop(paste("Attribute for absdiff() must be a numeric vector."))
  xm<-as.matrix(dist(x))
  cn<-paste('absdiff',attrname,sep=".")

  InitErgmm.latentcov(model,xm,cn,
                      mean=mean,var=var)
}
InitErgmm.nodematch<-function (model, attrname, diff=FALSE, mean=0, var=9) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
        stop(paste("nodematch() model term expected between 1 and 4 arguments, got ", 
                   nargs() - 1, sep = ""), call. = FALSE)

  x<-as.factor(model[["Yg"]] %v% attrname)
  if(!diff){
    xm<-as.matrix(dist(as.numeric(x)))==0
    cn<-paste('nodematch',attrname,sep=".")
    model<-InitErgmm.latentcov(model,xm,cn,mean=mean,var=var)
  }else{
    ls<-levels(x)
    xm<-as.matrix(dist(as.numeric(x)))==0
    mean<-rep(mean,length.out=length(ls))
    var<-rep(var,length.out=length(ls))
    for(li in 1:length(ls)){
      l<-ls[[li]]
      xm.l<-xm
      xm.l[x!=l,]<-FALSE
      xm.l[,x!=l]<-FALSE
      cn<-paste('nodematch',attrname,l,sep=".")
      model<-InitErgmm.latentcov(model,xm.l,cn,mean=mean[li],var=var[li])
    }
  }
  
  model
}
InitErgmm.sendercov<-function (model, attrname, force.factor=FALSE, mean=0, var=9) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("sendercov() model term expected between 1 and 4 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  if (!is.directed(model[["Yg"]]))
    stop("Sender covariates are not allowed with an undirected network; use 'socialitycov'", call.=FALSE)

  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname
  if(is.numeric(x) && ! force.factor){
    # Numeric covariate.
    xm<-matrix(x,n,n,byrow=FALSE)
    cn<-paste("sendercov",attrname,sep=".")
    model<-InitErgmm.latentcov(model,xm,cn,mean=mean,var=var)
  }else{
    # Factor covariate.
    x<-as.factor(x)
    ls<-levels(x)
    mean<-rep(mean,length.out=length(ls))
    var<-rep(var,length.out=length(ls))
    for(li in 2:length(ls)){
      l<-ls[[li]]
      xm<-matrix(x==l,n,n,byrow=FALSE)
      cn<-paste('sendercov',attrname,l,sep=".")
      model<-InitErgmm.latentcov(model,xm,cn,mean=mean[li],var=var[li])
    }
  }
  model
}
InitErgmm.receivercov<-function (model, attrname, force.factor=FALSE, mean=0, var=9) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("receivercov() model term expected between 1 and 3 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)
  if (!is.directed(model[["Yg"]]))
    stop("Receiver covariates are not allowed with an undirected network; use 'socialitycov'", call.=FALSE)

  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname
  if(is.numeric(x) && ! force.factor){
    # Numeric covariate.
    xm<-matrix(x,n,n,byrow=TRUE)
    cn<-paste("receivercov",attrname,sep=".")
    model<-InitErgmm.latentcov(model,xm,cn,mean=mean,var=var)
  }else{
    # Factor covariate.
    x<-as.factor(x)
    ls<-levels(x)
    mean<-rep(mean,length.out=length(ls))
    var<-rep(var,length.out=length(ls))
    for(li in 2:length(ls)){
      l<-ls[[li]]
      xm<-matrix(x==l,n,n,byrow=TRUE)
      cn<-paste('receivercov',attrname,l,sep=".")
      model<-InitErgmm.latentcov(model,xm,cn,mean=mean[li],var=var[li])
    }
  }
  model
}
InitErgmm.socialitycov<-function (model, attrname, force.factor=FALSE, mean=0, var=9) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("socialitycov() model term expected between 1 and 4 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  n<-network.size(model[["Yg"]])
  x<-model[["Yg"]] %v% attrname
  if(is.numeric(x) && ! force.factor){
    # Numeric covariate.
    xm<-matrix(x,n,n,byrow=FALSE)+matrix(x,n,n,byrow=TRUE)
    cn<-paste("socialitycov",attrname,sep=".")
    model<-InitErgmm.latentcov(model,xm,cn,mean=mean,var=var)
  }else{
    # Factor covariate.
    x<-as.factor(x)
    ls<-levels(x)
    mean<-rep(mean,length.out=length(ls))
    var<-rep(var,length.out=length(ls))
    for(li in 2:length(ls)){
      l<-ls[[li]]
      xm<-matrix(x==l,n,n,byrow=FALSE)+matrix(x==l,n,n,byrow=TRUE)
      cn<-paste('socialitycov',attrname,l,sep=".")
      model<-InitErgmm.latentcov(model,xm,cn,mean=mean[li],var=var[li])
    }
  }
  model
}
