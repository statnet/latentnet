InitErgmm.latent<-function(model, d, G=0, var.mul=1/30, var=NULL, var.df=3,
                           mean.var.mul=1/5, mean.var=NULL, pK=3){
  if (nargs()<2)
    stop(paste("latent() model term expected at least 1 argument, got ",
               nargs()-1, sep=""), call.=FALSE)
  model$d <- d
  model$G <- G

  model$prior$Z.var.mul<-var.mul
  model$prior$Z.var<-var
  model$prior$Z.var.df<-var.df
  model$prior$Z.mean.var.mul<-mean.var.mul
  model$prior$Z.mean.var<-mean.var
  model$prior$Z.pK<-pK
  
  model
}
