InitErgmm.latent<-function(model, d, G=0, var.mul=1/4, var=NULL, var.df.mul=1/2, var.df=NULL,
                           mean.var.mul=2, mean.var=NULL, pK.mul=1/2, pK=NULL){
  if (nargs()<2)
    stop(paste("latent() model term expected at least 1 argument, got ",
               nargs()-1, sep=""), call.=FALSE)
  model$d <- d
  model$G <- G

  model$prior$Z.var.mul<-var.mul
  model$prior$Z.var<-var
  model$prior$Z.var.df.mul<-var.df.mul
  model$prior$Z.var.df<-var.df
  model$prior$Z.mean.var.mul<-mean.var.mul
  model$prior$Z.mean.var<-mean.var
  model$prior$Z.pK.mul<-pK.mul
  model$prior$Z.pK<-pK
  
  model
}
