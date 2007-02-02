InitErgmm.absdiff<-function (model, attrname,
                             mean=0, var=9) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 2:4))
    stop(paste("absdiff() model term expected between 1 and 3 arguments, got ", 
                                   nargs() - 1, sep = ""), call. = FALSE)

  x<-model$Yg %v% attrname
  if(!is.vector(x)||!is.numeric(x))
    stop(paste("Attribute for absdiff() must be a numeric vector."))
  xm<-as.matrix(dist(x))
  cn<-paste('absdiff',attrname,sep=".")

  InitErgmm.latentcov(model,xm,cn,
                      mean=mean,var=var)
}
