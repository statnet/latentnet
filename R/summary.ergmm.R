summary.ergmm <- function (object, point.est=c("pmean","mkl"), quantiles=c(.025,.975),se=FALSE,...)
{
  extraneous.argcheck(...)
  ## Just for convenience.
  sample<-object$sample
  model<-object$model
  control<-object$control
  d<-model$d
  p<-model$p
  G<-model$G
  n<-network.size(model$Yg)
  sample.size<-control$sample.size
  

  summ<-list(ergmm=object,model=model)

  ## Compute the two-stage MLE point estimates.
  if("mle" %in% point.est){
    if(is.null(object$mle)){
      stop("MLE was not computed for this fit.")
    }
    if(is.null(object$mle$cov)){
      if(se){
        ## Refit the MLE (mostly for the Hessian)
        mle<-find.mle(model,object$mle,control=control,hessian=TRUE)
      }
      else mle<-object$mle

      if(d>0) mle$Z<-scale(mle$Z,scale=FALSE)
      if(p>0){
        if(se){
          beta.hess<-mle$hessian[1:p,1:p]
          
          beta.cov <- try(robust.inverse(-beta.hess), silent=TRUE)
          if(inherits(beta.cov,"try-error")){
            warning("Coefficient Hessian appears to be singular. Using a less accurate estimate.")
            beta.cov <- diag(1/diag(-beta.hess))
          }
        
        colnames(beta.cov) <- rownames(beta.cov) <- model$coef.names
        
        mle$cov<-beta.cov
        mle$cor<-cov2cor(beta.cov)
        
        z.values<-mle$beta/sqrt(diag(mle$cov))
        coef.table<-data.frame(mle$beta,sqrt(diag(mle$cov)),
                               z.values,pnorm(abs(z.values),0,1,lower.tail=FALSE),
                             row.names=model$coef.names)
        }
        else coef.table<-data.frame(mle$beta,
                                    row.names=model$coef.names)
        colnames(coef.table)<-c("Estimate",if(se) "Std. Error",if(se) "z value",if(se) "Pr(>|z|)")
        mle$coef.table<-coef.table
      }
      mle$llk.bic<- -2*mle$llk + (p+n*d)*log(sum(network.size(model$Yg)))
      if(G>0){
        mle<-c(mle,find.clusters(G,mle$Z))
        mle$mbc.bic<- -bic(if(d==1) "V" else "VII",mle$mbc.llk,n,d,G)
        mle$bic<-mle$llk.bic+mle$mbc.bic
      }
      else mle$bic<-mle$llk.bic

      object$mle<-mle
    }
    summ$mle<-object$mle
  }

  ## Compute the posterior mean point estimates.
  if("pmean" %in% point.est){
    if(is.null(object$sample)){
      stop("MCMC was not was not run for this fit.")
    }
    if(is.null(object$pmean)){
      pmean<-list()
      for(name in names(sample)){
        if(is.null(sample[[name]])) next
        name.dim<-length(dim(sample[[name]]))
        pmean[[name]]<-{
          if(name.dim<2) mean(sample[[name]])
          else if(name.dim==2) apply(sample[[name]],2,mean)
          else if(name.dim==3) apply(sample[[name]],2:3,mean)
        }
      }
      if(G>0){
        pmean$Z.pZK<-t(apply(sample$Z.K,2,tabulate,G))/sample.size
        pmean$Z.K<-apply(pmean$Z.pZK,1,which.max)
      }
    
      beta.cov<-cov(sample$beta)
      colnames(beta.cov) <- rownames(beta.cov) <- model$coef.names
      
      pmean$cov<-beta.cov
      pmean$cor<-cov2cor(beta.cov)
      
      beta.q0<-apply(sample$beta,2,function(x) mean(x<=0))
      
      coef.table<-data.frame(pmean$beta,
                           t(apply(sample$beta,2,function(x)quantile(x,quantiles))),
                           beta.q0,
                           row.names=model$coef.names)
      colnames(coef.table)<-c("Estimate",paste(quantiles*100,"%",sep=""),"Quantile of 0")
      pmean$coef.table<-coef.table
      object$pmean<-pmean
    }
    summ$pmean<-object$pmean
  }
  ## Compute the MKL point estimates.
  if("mkl" %in% point.est){
    if(is.null(object$mkl)){
      stop("MKL was not produced for this fit.")
    }
    mkl<-summ$mkl<-object$mkl
    coef.table<-data.frame(mkl$beta,
                           row.names=model$coef.names)
    colnames(coef.table)<-c("Estimate")
    summ$mkl$coef.table<-coef.table
  }

  if("pmode" %in% point.est){
    if(is.null(object$pmode)){
      stop("Conditional posterior mode was not computed for this fit.")
    }
    summ$pmode<-object$pmode
  }

  class(summ)<-'summary.ergmm'
  summ
}

print.summary.ergmm<-function(x,...){
  ## For convenience
  model<-x$model
  control<-x$ergmm$control
  
  cat("\n==========================\n")
  cat("Summary of model fit\n")
  cat("==========================\n\n")

  cat("Formula:   ")
  print(model$formula)
  cat("Attribute: ")
  if(is.null(model$response)) cat("edges") else cat(model$response)
  cat("\n")
  cat("Model:    ",model$family,"\n")
  
  digits = max(3, getOption("digits") - 3)
  
  cat ("MCMC sample of size ", control$sample.size, ", draws are ",
       control$interval," iterations apart, after burnin of ",control$burnin, " iterations.\n",sep="")
       
  if(!is.null(x$pmean)){
    cat("Covariate coefficients posterior means:\n")
    print(x$pmean$coef.table)
    cat("\n")
  }

  if(!is.null(x$mle)){
    cat("Covariate coefficients MLE:\n")
    print(x$mle$coef.table)
    cat("\n")

    cat("Overall BIC:       ", x$mle$bic,"\n")
    if(model$G>0){
      cat("Likelihood BIC:    ", x$mle$llk.bic,"\n")
      cat("Clustering BIC:    ", x$mle$mbc.bic,"\n")
    }
    cat("\n")
  }

  if(!is.null(x$mkl)){
    cat("Covariate coefficients MKL:\n")
    print(x$mkl$coef.table)
    cat("\n\n")
  }
  if(!is.null(x$pmode)){
    cat("Covariate coefficients posterior mode:\n")
    print(x$pmode$coef.table)
    cat("\n\n")
  }
}
