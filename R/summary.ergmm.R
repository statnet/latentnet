summary.ergmm <- function (ergmm.fit, point.est=c("pmean","mkl"), quantiles=c(.025,.975),se=TRUE)
{
  ## Just for convenience.
  samples<-ergmm.fit$samples
  model<-ergmm.fit$model
  control<-ergmm.fit$control
  d<-model$d
  p<-model$p
  G<-model$G
  n<-network.size(model$Yg)
  samplesize<-control$samplesize
  

  summ<-list(ergmm=ergmm.fit,model=model)

  ## Compute the two-stage MLE point estimates.
  if("mle" %in% point.est){
    if(is.null(ergmm.fit$mle$cov)){
      if(model$sender || model$receiver || model$sociality)
        warning("Fitting random effects as fixed effects.")
      if(se){
        ## Refit the MLE (mostly for the Hessian)
        mle<-find.mle(model,ergmm.fit$mle,control=control,hessian=TRUE)
      }
      else mle<-ergmm.fit$mle

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
      mle$llk.bic<- -2*mle$llk + (p+n*d + (model$sender + model$receiver + model$sociality)*n )*log(sum(network.size(model$Yg)))
      if(G>0){
        mle<-c(mle,find.clusters(G,mle$Z))
        mle$mbc.bic<- -bic(if(d==1) "V" else "VII",mle$mbc.llk,n,d,G)
        mle$bic<-mle$llk.bic+mle$mbc.bic
      }
      else mle$bic<-mle$llk.bic

      if(model$sender){
        if(model$intercept){
          mle$beta[1]<-mle$beta[1]-mean(mle$sender)
          mle$sender<-mle$sender-mean(mle$sender)
        }
        mle$sender.var<-mean(mle$sender^2)
      }
        
      if(model$receiver){
        if(model$intercept){
          mle$beta[1]<-mle$beta[1]-mean(mle$receiver)
          mle$receiver<-mle$receiver-mean(mle$receiver)
        }
        mle$receiver.var<-mean(mle$receiver^2)
      }

      if(model$sociality){
        if(model$intercept){
          mle$beta[1]<-mle$beta[1]-mean(mle$sociality)
          mle$sociality<-mle$sociality-mean(mle$sociality)
        }
        mle$sociality.var<-mean(mle$sociality^2)
      }

      ergmm.fit$mle<-mle
    }
    summ$mle<-ergmm.fit$mle
  }

  ## Compute the posterior mean point estimates.
  if("pmean" %in% point.est){
    if(is.null(ergmm.fit$pmean)){
      pmean<-list()
      for(name in names(samples)){
        if(is.null(samples[[name]])) next
        name.dim<-length(dim(samples[[name]]))
        pmean[[name]]<-{
          if(name.dim<2) mean(samples[[name]])
          else if(name.dim==2) apply(samples[[name]],2,mean)
          else if(name.dim==3) apply(samples[[name]],2:3,mean)
        }
      }
      if(G>0){
        pmean$Z.pZK<-t(apply(samples$Z.K,2,tabulate,G))/samplesize
        pmean$Z.K<-apply(pmean$Z.pZK,1,which.max)
      }
    
      beta.cov<-cov(samples$beta)
      colnames(beta.cov) <- rownames(beta.cov) <- model$coef.names
      
      pmean$cov<-beta.cov
      pmean$cor<-cov2cor(beta.cov)
      
      beta.q0<-apply(samples$beta,2,function(x) mean(x<=0))
      
      coef.table<-data.frame(pmean$beta,
                           t(apply(samples$beta,2,function(x)quantile(x,quantiles))),
                           beta.q0,
                           row.names=model$coef.names)
      colnames(coef.table)<-c("Estimate",paste(quantiles*100,"%",sep=""),"Quantile of 0")
      pmean$coef.table<-coef.table
      ergmm.fit$pmean<-pmean
    }
    summ$pmean<-ergmm.fit$pmean
  }
  ## Compute the MKL point estimates.
  if("mkl" %in% point.est){
    mkl<-summ$mkl<-ergmm.fit$mkl
    coef.table<-data.frame(mkl$beta,
                           row.names=model$coef.names)
    colnames(coef.table)<-c("Estimate")
    summ$mkl$coef.table<-coef.table
  }

  if("pmode" %in% point.est){
    summ$pmode<-ergmm.fit$pmode
  }

  class(summ)<-'summary.ergmm'
  summ
}

print.summary.ergmm<-function(summ,...){
  ## For convenience
  model<-summ$model
  control<-summ$control
  
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
  
  cat ("MCMC sample of size ", control$samplesize, ", samples are ",
       control$interval," iterations apart, after burnin of ",control$burnin, " iterations.\n",sep="")
       
  if(!is.null(summ$pmean)){
    cat("Covariate coefficients posterior means:\n")
    print(summ$pmean$coef.table)
    cat("\n")
    if(!is.null(summ$pmean$sender.var))
      cat("Sender effect variance: ",summ$pmean$sender.var,".\n", sep="")
    if(!is.null(summ$pmean$receiver.var))
      cat("Receiver effect variance: ",summ$pmean$receiver.var,".\n", sep="")
    if(!is.null(summ$pmean$sociality.var))
      cat("Sociality effect variance: ",summ$pmean$sociality.var,".\n", sep="")
  }

  if(!is.null(summ$mle)){
    cat("Covariate coefficients MLE:\n")
    print(summ$mle$coef.table)
    cat("\n")

    cat("Overall BIC:       ", summ$mle$bic,"\n")
    if(model$G>0){
      cat("Likelihood BIC:    ", summ$mle$llk.bic,"\n")
      cat("Clustering BIC:    ", summ$mle$mbc.bic,"\n")
    }
    cat("\n")
  }
  if(!is.null(summ$pmean)){
    cat("Covariate coefficients posterior mean:\n")
    print(summ$pmean$coef.table)
    cat("\n\n")
  }
  if(!is.null(summ$mkl)){
    cat("Covariate coefficients MKL:\n")
    print(summ$mkl$coef.table)
    cat("\n\n")
  }
  if(!is.null(summ$pmode)){
    cat("Covariate coefficients posterior mode:\n")
    print(summ$pmode$coef.table)
    cat("\n\n")
  }
}
