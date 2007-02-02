InitErgmm.rsender<-function(model, var=1, var.df=3){
  if (!is.directed(model$Yg))
    stop("Sender effects are not allowed with an undirected network; use 'sociality'", call.=FALSE)
  model$sender<-TRUE
  model$prior$sender.var<-var
  model$prior$sender.var.df<-var.df
  model
}
