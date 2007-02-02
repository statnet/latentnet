InitErgmm.socialitycov<-function (model, attrname, force.factor=FALSE, mean=0, var=9) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
    stop(paste("socialitycov() model term expected between 1 and 4 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  n<-network.size(model$Yg)
  x<-model$Yg %v% attrname
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
