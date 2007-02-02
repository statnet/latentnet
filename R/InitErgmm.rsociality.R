InitErgmm.rsociality<-function(model, var=1, var.df=3){
  model$sociality<-TRUE
  model$prior$sociality.var<-var
  model$prior$sociality.var.df<-var.df
  model
}
