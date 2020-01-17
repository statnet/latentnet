#  File R/print.ergmm.model.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
#' @export
print.ergmm.model<-function(x,...){
  cat("Exponential Random Graph Mixed Model definition\n")
  cat("Formula: ")
  print(x[["formula"]])
  if(!is.null(x[["response"]]) && mode(x[["response"]])=="character")
    cat("Attribute of interest:",x[["response"]],"\n")
  cat("Family:",x[["family"]],"\n")

  cat("Terms:\n")
  if(length(x[["coef.names"]])){
    cat("- fixed effects:\n")
    cat(paste(" - ",x[["coef.names"]],sep="",collapse="\n"),"\n")
  }
  if(x[["d"]]>0){
    cat("- latent space of",x[["d"]],"dimensions")
    if(x[["G"]]>0){
      cat(" in",x[["G"]],"clusters\n")
    }else cat("\n")
  }
  if(x[["sender"]]) cat("- random sender effects\n")
  if(x[["receiver"]]) cat("- random receiver effects\n")
  if(x[["sociality"]]) cat("- random sociality effects\n")
  if(x[["dispersion"]]) cat("- dispersion parameter\n")
}
