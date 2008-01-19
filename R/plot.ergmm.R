"plot.ergmm" <- function (x, ..., vertex.cex=1, vertex.sides=16*ceiling(sqrt(vertex.cex)),
                          what="mkl",
                          main = NULL, xlab=NULL, ylab=NULL, xlim=NULL,ylim=NULL,
                          object.scale=formals(plot.network.default)$object.scale,
                          pad=formals(plot.network.default)$pad,
                          cluster.col=c("red","green","blue","cyan","magenta","orange","yellow","purple"),
                          vertex.col=NULL, print.formula=TRUE,
                          edge.col=8,
                          pie = FALSE,
                          labels=FALSE,
                          plot.means=TRUE,plot.vars=TRUE,
                          suppress.axes=FALSE,
                          jitter1D=1,curve1D=TRUE,suppress.center=FALSE,
                          density.par=list())
{

  ## For convenience...
  Yg<-x$model$Yg
  distances<-NULL
  n<-network.size(Yg)
  d<-x$model$d
  G<-x$model$G
  if(G<1) pie<-FALSE

  ## Set default axis labels.
  if(d==1){
    if (is.null(xlab)) 
      xlab <- ""
    if (is.null(ylab)) 
      ylab <- ""
  }else if(d==2){    
    if (is.null(xlab)) 
      xlab <- expression(Z[1])
    if (is.null(ylab)) 
      ylab <- expression(Z[2])
  }else if(d>2){
    if (is.null(xlab)) 
      xlab <- "First principal component of Z"
    if (is.null(ylab)) 
      ylab <- "Second principal component of Z"
  }

  ## Find the requested plotting coordinates.
  ## Some "requests" require a substantially different code path, unfortunately.
  if(class(what)=="ergmm.par"){
    summ<-what
    Z.pos <- summ$Z
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    Z.pZK<-summ$Z.pZK
    if (is.null(main)) 
      main <- paste(deparse(substitute(what))," Latent Positions of ", 
                    deparse(substitute(x)),sep="")

  }else if(what=="start" || what=="burnin.start"){
    summ<-x$start
    Z.pos <- summ$Z
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (is.null(main)) 
      main <- paste("Initial Latent Positions of ", 
                    deparse(substitute(x)),sep="")

  }else if(what=="sampling.start"){
    summ<-x$sampling.start
    Z.pos <- summ$Z
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (is.null(main)) 
      main <- paste("Post-Burn-In Latent Positions of ", 
                    deparse(substitute(x)),sep="")

  }else if(what=="mle"){
    summ<-summary(x,point.est=c("mle"),se=FALSE)
    Z.pos <- summ$mle$Z
    summ<-summ$mle
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    plot.means<-plot.vars<-FALSE
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (is.null(main)) 
      main <- paste("Multistage MLEs of Latent Positions of", 
                    deparse(substitute(x)))

  }else if(what=="pmean"){
    summ<-summary(x,point.est=c("pmean"))
    Z.pos <- summ$pmean$Z
    summ<-summ$pmean
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    Z.pZK<-summ$Z.pZK
    if (is.null(main)) 
      main <- paste("Posterior Mean Positions of", 
                    deparse(substitute(x)))

  }else if(what=="mkl"){
    summ<-summary(x,point.est=c("pmean","mkl"))
    Z.pos <- summ$mkl$Z
    if(!is.null(x$mkl$mbc)){
      Z.mean<-summ$mkl$mbc$Z.mean
      Z.var<-summ$mkl$mbc$Z.var
    }else{
      if(!is.null(summ$pmean$Z.mean)) Z.mean<-summ$pmean$Z.mean
      else plot.means<-FALSE
      Z.var<-summ$pmean$Z.var
    }
    Z.K<-summ$pmean$Z.K
    Z.pZK<-summ$pmean$Z.pZK
    summ<-summ$mkl
    if (is.null(main)) 
      main <- paste("MKL Latent Positions of", 
                    deparse(substitute(x)))

  }else if(what=="pmode"){
    summ<-summary(x,point.est=c("pmode"))
    Z.pos <- summ$pmode$Z
    summ<-summ$pmode
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (is.null(main)) 
      main <- paste("Posterior Mode Latent Positions of", 
                    deparse(substitute(x)))

  }else if(what=="cloud"){
    summ<-summary(x,point.est=c("pmean","mkl"))
    Z.pos <- summ$mkl$Z
    if(d!=2) stop("Cloud plots are only available for 2D latent space models.")
    
    if(!is.null(x$mkl$mbc)){
      Z.mean<-summ$mkl$mbc$Z.mean
      Z.var<-summ$mkl$mbc$Z.var
    }else{
      if(!is.null(summ$pmean$Z.mean)) Z.mean<-summ$pmean$Z.mean
      else plot.means<-FALSE
      Z.var<-summ$pmean$Z.var
    }
    Z.K<-summ$pmean$Z.K
    Z.pZK<-summ$pmean$Z.pZK
    summ<-summ$mkl
    if (is.null(main)) 
      main <- paste("MKL Latent Positions of", 
                    deparse(substitute(x)))
    plot(matrix(c(x$samples$Z),ncol=2),pch=".")
    points(Z.pos,col=cluster.col[Z.K])
    points(Z.mean,col=cluster.col)
    
    return(invisible(NULL))

  }else if(what=="density"){
    if(is.null(density.par$totaldens)) density.par$totaldens <- TRUE
    if(is.null(density.par$subdens)) density.par$subdens <- TRUE
    if(is.null(density.par$mfrow)){
      wanted<-density.par$totaldens+density.par$subdens*G
      density.par$mfrow<-rep(min(ceiling(sqrt(wanted)),4),2)
    }
    
    summ<-summary(x,point.est=c("pmean","mkl"))
    Z.pos <- summ$mkl$Z
    if(d!=2) stop("Density plots are only available for 2D latent space models.")

    if(!require(KernSmooth,quietly=TRUE)){
      stop("The 'density' option requires the 'KernSmooth' package.")
    }

    old.par<-par(mfrow=density.par$mfrow,mar=c(2.5,2.5,1,1))

    Z.all<-matrix(c(aperm(x$samples$Z,c(2,1,3))),ncol=2)

    if(density.par$totaldens){
      plot(Z.all,type='n',xlab=xlab,ylab=ylab,...)
      title(main=paste("Posterior density of",deparse(substitute(x))), cex.main=0.7, ...)
      Z.bkde <- bkde2D(Z.all,0.2,c(201,201))
      image(Z.bkde$x1,Z.bkde$x2,Z.bkde$fhat,col=grey(seq(1,0,length=255)),add=TRUE)
      box()
    }
    
    if(G>1 && density.par$subdens){
      Z.K.all <- c(t(x$samples$Z.K))
      for(i in 1:G){
        plot(Z.all,main=paste("Class",i),type="n",...)
        Z.bkde <- bkde2D(Z.all[Z.K.all==i,],0.2,c(101,101))
        col<-c(col2rgb(cluster.col[i])/255)
        image(Z.bkde$x1,Z.bkde$x2,Z.bkde$fhat,add=TRUE,
              col=rgb(seq(1,col[1],length=255),
                seq(1,col[2],length=255),
                seq(1,col[3],length=255)))
        contour(Z.bkde$x1,Z.bkde$x2,Z.bkde$fhat,add=TRUE, nlevels=4,
                drawlabels=FALSE,
                col="white")
        box()
      }
    }
    
    par(old.par)

    return(invisible(NULL))
    
  }else if(is.numeric(what) && round(what)==what){
    summ<-x$samples[[what]]
    Z.pos <- summ$Z
    Z.mean<-summ$Z.mean
    Z.var<-summ$Z.var
    Z.K<-summ$Z.K
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (is.null(main)) 
      main <- paste("Iteration #",what," Latent Positions of ", 
                    deparse(substitute(x)),sep="")
  }else stop("Invalid latent space position estimate type.")

  ## Transform coordinates for dimensions other than 2.
  if(d==1){    
    Z.pos<-coords.1D(Z.pos,curve1D,jitter1D)
    if(curve1D){
      distances<-as.matrix(dist(Z.pos))
      distances<-distances/max(distances)
    }
  } else if(d>2){ ## I.e. high latent dimension.
    ## Plot the first two principal components.
    prc<-prcomp(Z.pos)
    Z.pos<-predict(prc,Z.pos)[,1:2]
  }

  ## Set default vertex color.
  if(is.null(vertex.col)){
    if(is.latent.cluster(x) && !pie)
      vertex.col <- cluster.col[Z.K]
    else vertex.col<-cluster.col[1]
  }
  else if(length(vertex.col)==1 && is.character(vertex.col)){
    trycol <- as.numeric(as.factor(unlist(get.vertex.attribute(Yg,vertex.col))))
    if(!all(is.na(trycol))){
      vertex.col <- cluster.col[trycol]
    }
  }

  ## Find the bounds of the plotting region.
  xylim<-ergmm.plotting.region(Z.pos,if(plot.means) Z.mean,if(plot.vars) Z.var,!suppress.center,pad)
  if(is.null(xlim)) xlim<-xylim$xlim else xylim$xlim<-xlim
  if(is.null(ylim)) ylim<-xylim$ylim else xylim$ylim<-ylim
    
  ## Go forth, and plot the network.
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
               displaylabels=labels&&!(x$model$d==1 && curve1D==TRUE),
               tick=!(x$model$d==1 && curve1D==TRUE),
               edge.col=edge.col,
               ...)
  
  options(warn=old.warn)
  
  ## For 1D plots, only plot horizontal axis ticks.
  if(x$model$d==1 && curve1D==TRUE) {
    axis(1)
  }

  ## Add the model formula as a subtitle.
  if(print.formula){
    fmla <- x$model$formula   
    xformula <- paste(fmla[2],fmla[1],fmla[-c(1:2)],collapse=" ")
    if(!is.null(x$model$response)) xformula<-paste(xformula,"   (",x$model$response,")",sep="")
    title(main = xformula, line = 1, cex.main = 0.7)
  }

  ## Plot pie charts.
  if(pie){
    piesize<-rep(ergmm.plotting.vertex.radius(vertex.cex,xylim,object.scale),length=n)
    pie.order<-order(piesize,decreasing=TRUE)
    for(i in 1:n){
      ergmm.drawpie(Z.pos[pie.order[i],],piesize[pie.order[i]],Z.pZK[pie.order[i],],n=50,cols=cluster.col)
    }
  }

  ## Mark the center.
  if(!suppress.center)
    points(cbind(0,0),pch="+")
  Z.mean<-if(G>0)Z.mean else cbind(0,0)
  if(x$model$d==1)
    Z.mean<-coords.1D(Z.mean,curve1D,jitter1D)
  else if(x$model$d>2){
    Z.mean<-predict(prc,Z.mean)
  }

  ## Mark the cluster means.
  if(plot.means)
    points(Z.mean,pch="+",col=cluster.col)

  ## Plot the cluster standard deviations.
  if(plot.vars)
      symbols(Z.mean,circles=sqrt(Z.var),fg=cluster.col,inches=FALSE,add=TRUE,asp=1)
  
  invisible(Z.pos)
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
