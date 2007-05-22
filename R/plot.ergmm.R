"plot.ergmm" <- function (ergmm.fit, ..., vertex.cex=1, vertex.sides=16*ceiling(sqrt(vertex.cex)),
                          which.par="mkl",
                          main = NULL, xlab=NULL, ylab=NULL, xlim=NULL,ylim=NULL,
                          object.scale=formals(plot.network.default)$object.scale,
                          pad=formals(plot.network.default)$pad,
                          cluster.col=c("red","green","blue","cyan","magenta","orange","yellow","purple"),
                          vertex.col=NULL, print.formula=TRUE,
                          edge.col=8,
                          pie = FALSE,
                          labels=TRUE,
                          rand.eff=NULL,
                          plot.means=TRUE,plot.vars=TRUE,
                          suppress.axes=FALSE,
                          jitter1D=1,curve1D=TRUE,suppress.center=FALSE)
{

  Yg<-ergmm.fit$model$Yg
  Ym <- sociomatrix(Yg)
  distances<-NULL

  n<-network.size(Yg)

  if (missing(xlab)) 
    xlab <- ""
  if (missing(ylab)) 
    ylab <- ""
  
  if(class(which.par)=="ergmm.par"){
    summ<-which.par
    Z.pos <- summ$Z
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    Z.pZK<-summ$Z.pZK
    if (missing(main)) 
      main <- paste(deparse(substitute(which.par))," Latent Positions of ", 
                    deparse(substitute(ergmm.fit)),sep="")
  }else if(which.par=="start"){
    summ<-ergmm.fit$start
    Z.pos <- summ$Z
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (missing(main)) 
      main <- paste("Initial Latent Positions of ", 
                    deparse(substitute(ergmm.fit)),sep="")
  }else if(which.par=="mle"){
    summ<-summary(ergmm.fit,point.est=c("mle"),se=FALSE)
    Z.pos <- summ$mle$Z
    summ<-summ$mle
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    plot.means<-plot.vars<-FALSE
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (missing(main)) 
      main <- paste("Multistage MLEs of Latent Positions of", 
                    deparse(substitute(ergmm.fit)))
  }else if(which.par=="pmean"){
    summ<-summary(ergmm.fit,point.est=c("pmean"))
    Z.pos <- summ$pmean$Z
    summ<-summ$pmean
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    Z.pZK<-summ$Z.pZK
    if (missing(main)) 
      main <- paste("Posterior Mean Positions of", 
                    deparse(substitute(ergmm.fit)))
  }else if(which.par=="mkl"){
    summ<-summary(ergmm.fit,point.est=c("pmean","mkl"))
    Z.pos <- summ$mkl$Z
    Z.mean<-summ$mkl$mbc$Z.mean
    Z.var<-summ$mkl$mbc$Z.var
    Z.K<-summ$pmean$Z.K
    Z.pZK<-summ$pmean$Z.pZK
    summ<-summ$mkl
    if (missing(main)) 
      main <- paste("MKL Latent Positions of", 
                    deparse(substitute(ergmm.fit)))
  }else if(which.par=="pmode"){
    summ<-summary(ergmm.fit,point.est=c("pmode"))
    Z.pos <- summ$pmode$Z
    summ<-summ$pmode
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (missing(main)) 
      main <- paste("Posterior Mode Latent Positions of", 
                    deparse(substitute(ergmm.fit)))
  }else if(which.par=="cloud"){
    summ<-summary(ergmm.fit,point.est=c("pmean","mkl"))
    Z.pos <- summ$mkl$Z
    Z.mean<-summ$mkl$mbc$Z.mean
    Z.var<-summ$mkl$mbc$Z.var
    Z.K<-summ$pmean$Z.K
    Z.pZK<-summ$pmean$Z.pZK
    summ<-summ$mkl
    if (missing(main)) 
      main <- paste("MKL Latent Positions of", 
                    deparse(substitute(ergmm.fit)))
    plot(matrix(c(ergmm.fit$samples$Z),ncol=2),pch=".")
    points(Z.pos,col=cluster.col[Z.K])
    points(Z.mean,col=cluster.col)
    invisible(NULL)
  }else if(is.numeric(which.par) && round(which.par)==which.par){
    summ<-ergmm.fit$samples[[which.par]]
    Z.pos <- summ$Z
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (missing(main)) 
      main <- paste("Iteration #",which.par," Latent Positions of ", 
                    deparse(substitute(ergmm.fit)),sep="")
  }else stop("Invalid latent space position estimate type.")
  
  if(dim(Z.pos)[2]==1){
    Z.pos<-coords.1D(Z.pos,curve1D,jitter1D)
    if(curve1D){
      distances<-as.matrix(dist(Z.pos))
      distances<-distances/max(distances)
    }
  }
  else
    if(dim(Z.pos)[2]>2){ ## I.e. high latent dimension.
      ## Plot the first two principal components.
      prc<-prcomp(Z.pos)
      Z.pos<-predict(prc,Z.pos)[,1:2]
    }
  
  if(is.null(vertex.col)){
    if(is.latent.cluster(ergmm.fit) && !pie)
      vertex.col <- cluster.col[Z.K]
    else vertex.col<-cluster.col[1]
  }
  else if(length(vertex.col)==1 && is.character(vertex.col)){
    trycol <- as.numeric(as.factor(unlist(get.vertex.attribute(Yg,vertex.col))))
    if(!all(is.na(trycol))){
      vertex.col <- cluster.col[trycol]
    }
  }
  
  if(!missing(rand.eff) && (rand.eff[1]=="total" || ergmm.fit$model[rand.eff[1]][[1]])){
    if(rand.eff=="total")
      rand.eff.mul<-exp((summ["sender"][[1]]+summ["receiver"][[1]])/2)
    else      
      rand.eff.mul<-exp(summ[rand.eff][[1]]/2)
    rand.eff.mul<-rand.eff.mul/mean(rand.eff.mul)
    vertex.cex<-vertex.cex*rand.eff.mul
  }
  
  xylim<-ergmm.plotting.region(Z.pos,if(plot.means) Z.mean,if(plot.vars) Z.var,!suppress.center,pad)
  if(is.null(xlim)) xlim<-xylim$xlim else xylim$xlim<-xlim
  if(is.null(ylim)) ylim<-xylim$ylim else xylim$ylim<-ylim
    
  old.warn<-options()$warn
  options(warn=-1)
  
  plot.network(Yg,coord=Z.pos,
               main=main,
               vertex.col=if(is.null(vertex.col)) 0 else vertex.col,
               vertex.cex=vertex.cex,
               xlab=xlab,
               ylab=ylab,
               xlim=xlim,
               ylim=ylim,
               object.scale=object.scale,
               pad=pad,
               suppress.axes=suppress.axes,
               vertex.sides=16*sqrt(vertex.cex),
               jitter=FALSE,
               usecurve=curve1D,
               edge.curve=distances,
               displaylabels=labels&&!(ergmm.fit$model$d==1 && curve1D==TRUE),
               tick=!(ergmm.fit$model$d==1 && curve1D==TRUE),
               edge.col=edge.col,
               ...)
  options(warn=old.warn)
  
  ## For 1D plots, only plot horizontal axis ticks.
  if(ergmm.fit$model$d==1 && curve1D==TRUE) {
    axis(1)
  }
  
  if(print.formula){
    aaa <- ergmm.fit$model$formula   
    xformula <- paste(aaa[2],aaa[1],aaa[-c(1:2)],collapse=" ")
    title(main = xformula, line = 1, cex.main = 0.7)
  }
  
  if(pie){
    piesize<-rep(ergmm.plotting.vertex.radius(vertex.cex,xylim,object.scale),length=n)
    for(i in 1:n){
      ergmm.drawpie(Z.pos[i,],piesize[i],Z.pZK[i,],n=50,cols=cluster.col)
    }
  }
  if(!suppress.center)
    points(cbind(0,0),pch="+")
  Z.mean<-if(ergmm.fit$model$G>0)Z.mean else cbind(0,0)
  if(ergmm.fit$model$d==1)
    Z.mean<-coords.1D(Z.mean,curve1D,jitter1D)
  else if(ergmm.fit$model$d>2){
    Z.mean<-predict(prc,Z.mean)
  }
  if(plot.means)
    points(Z.mean,pch="+",col=cluster.col)
  if(plot.vars)
      symbols(Z.mean,circles=sqrt(Z.var),fg=cluster.col,inches=FALSE,add=TRUE,asp=1)
  
  invisible(NULL)
}

coords.1D<-function(Z,curve1D,jitter1D){
  if(curve1D){
    Z.pos<-cbind(Z,rep(0,length(Z)))
  }else{
    jitter<-rnorm(length(Z),0,jitter1D*sd(Z)/sqrt(length(Z)))
    Z.pos<-cbind(Z,Z)+cbind(jitter,-jitter)
  }
  Z.pos
}

ergmm.plotting.region<-function(Z.pos,Z.mean,Z.var,include.center,pad){
  xylim<-t(apply(rbind(Z.pos+pad,Z.pos-pad,
                       Z.mean+pad,Z.mean-pad,
                       if(!is.null(Z.mean) && !is.null(Z.var)) Z.mean+sqrt(Z.var)+pad,
                       if(!is.null(Z.mean) && !is.null(Z.var)) Z.mean-sqrt(Z.var)-pad,
                       if(include.center) 0+pad,
                       if(include.center) 0-pad,
                       if(is.null(Z.mean) && length(Z.var)==1) matrix(sqrt(Z.var)+pad,ncol=2,byrow=FALSE),
                       if(is.null(Z.mean) && length(Z.var)==1) matrix(-sqrt(Z.var)-pad,ncol=2,byrow=FALSE)
                       ),2,range))
  xlim<-xylim[1,]
  ylim<-xylim[2,]

  ### Shamelessly stolen from Carter's plot.network.default code.
  xrng<-diff(xlim)          #Force scale to be symmetric
  yrng<-diff(ylim)
  xctr<-(xlim[2]+xlim[1])/2 #Get center of plotting region
  yctr<-(ylim[2]+ylim[1])/2
  if(xrng<yrng)
    xlim<-c(xctr-yrng/2,xctr+yrng/2)
  else
    ylim<-c(yctr-xrng/2,yctr+xrng/2)
  ### End stolen part
  
  list(xlim=xlim,ylim=ylim)
}

ergmm.plotting.vertex.radius<-function(vertex.cex,xylim,object.scale){
  baserad<-min(diff(xylim$xlim),diff(xylim$ylim))*object.scale
  vertex.cex*baserad
}
