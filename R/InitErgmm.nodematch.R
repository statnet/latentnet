InitErgmm.nodematch<-function (model, attrname, diff=FALSE, mean=0, var=9) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:5))
        stop(paste("nodematch() model term expected between 1 and 4 arguments, got ", 
                   nargs() - 1, sep = ""), call. = FALSE)

  x<-as.factor(model$Yg %v% attrname)
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
