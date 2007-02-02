InitErgmm.rreceiver<-function(model, var=1, var.df=3){
  if (!is.directed(model$Yg))
    stop("receiver effects are not allowed with an undirected network; use 'sociality'", call.=FALSE)
  model$receiver<-TRUE
  model$prior$receiver.var<-var
  model$prior$receiver.var.df<-var.df
  model
}
