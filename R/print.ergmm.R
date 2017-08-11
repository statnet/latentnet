#' @export
print.ergmm<-function(x,...){
cat("Fitted ERGMM model:\n")
print(x$model)
}

#' @export
show.ergmm <- print.ergmm
