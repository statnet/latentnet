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
    xm<-as.matrix.network(g,matrix.type="adjacency",x)
    cn<-x
  }else{
    xm<-as.matrix(x)
    cn<-if(!is.null(attrname)) attrname else paste("matrix",length(model$X)+1)
  }

  model$p<-model$p+1
  model$X<-c(model$X,list(xm))
  model$coef.names<-c(model$coef.names, cn)

  model$prior$beta.mean<-c(model$prior$beta.mean,mean)
  model$prior$beta.var<-c(model$prior$beta.var,var)

  model
}
