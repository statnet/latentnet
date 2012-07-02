#
# Return TRUE iff object x is a ergmm fit object
# or a latent model
#
is.latent<-function(x) inherits(x,"ergmm") && x[["model"]][["d"]]>0

#
# Return TRUE iff object x is a ergmm fit object
# or a latent model and a "latentcluster" model or fit
#
is.latent.cluster<-function(x) inherits(x,"ergmm") && x[["model"]][["d"]]>0 && x[["model"]][["G"]]>0

#
# Return TRUE iff object x is a ergmm fit object
# or a latent model and a sender random effect
#
is.sender<-function(x) inherits(x,"ergmm") && x[["model"]][["sender"]]
is.receiver<-function(x) inherits(x,"ergmm") && x[["model"]][["receiver"]]

extraneous.argcheck<-function(...){
  ## Because #$^$% R wants various functions implementing generics to
  ## functions to take ..., which is wonderful for missing spelling
  ## errors.
  
  if(length(list(...)))stop("Extraneous arguments passed: ",
                            paste(list(...)))
}


thin.ergmm<-function(x,by,...){
  extraneous.argcheck(...)
  if(x[["control"]][["threads"]]>1) warning("Multithreaded run output. Stuff might be broken.")
  S<-x[["control"]][["sample.size"]]
  s.kept<-seq(from=1,to=S,by=by)
  x[["sample"]]<-x[["sample"]][s.kept]
  x[["control"]][["interval"]]<-x[["control"]][["interval"]]*by
  x[["control"]][["sample.size"]]<-length(s.kept)
  x
}

xtabs.ergmm<-function(x,ref,min.plurality=0){
  ref->Reference
  apply(attr(x[["sample"]],"Q"),1,which.max)->Fitted
  admit<-sapply(seq_along(Fitted), function(i) attr(x[["sample"]],"Q")[i,Fitted[i]]>=min.plurality)
  Reference<-Reference[admit]
  Fitted<-Fitted[admit]
  xtabs(~Reference+Fitted)
}

clust.homogeneity<-function(x,ref,soft=TRUE,marg=FALSE){
  if(soft){
    G<-x[["model"]][["G"]]
    p.K<-sapply(sort(unique(ref)),function(g.ref){
      Z.K.g.ref<-x[["sample"]][["Z.K"]][,ref==g.ref,drop=FALSE]
      n.g.ref<-dim(Z.K.g.ref)[2]
      mean(sapply(2:n.g.ref,function(i) mean(sapply(1:(i-1),function(j) mean(Z.K.g.ref[,i]==Z.K.g.ref[,j])))))
    }
           )
  }else{
    p.K<-apply(xtabs.ergmm(x,ref),1,function(y)sum((y/sum(y))^2))    
  }

  if(marg) weighted.mean(p.K,tabulate(ref)*(tabulate(ref)-1))
  else p.K
}

## Concatenate a character list with commas and ands in the right places.
.paste.and <- function(x, oq='', cq=''){
  x <- paste(oq, x, cq, sep='')
  if(length(x)==0) return('')
  if(length(x)==1) return(x)
  if(length(x)==2) return(paste(x[1],'and',x[2]))
  if(length(x)>=3) return(paste(paste(x[-length(x)], collapse=", "),', and ',x[length(x)],sep=''))
}

## If EXPR is NULL, return NULLV, otherwise return EXPR.
NVL <- function(EXPR, NULLV) if(!is.null(EXPR)) EXPR else NULLV
  
