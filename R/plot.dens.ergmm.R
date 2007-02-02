"plot.dens.ergmm" <- function (ergmm.fit, ..., mle=FALSE, comp.mat = NULL,
                               label = NULL, label.col = "black",
                               xlab, ylab, main, label.cex = 0.8, edge.lwd = 1,
                               edge.col="gray", al = 0.1,
                               contours=0, density=FALSE, only.subdens = FALSE, 
                               drawarrows=FALSE,print.formula=TRUE,
                               contour.color=1, plotgraph=FALSE, pie = FALSE, piesize=0.07,
                               vertex.col=0,vertex.pch=19,vertex.cex=2,
                               cluster.col=c("red","green","blue","cyan",
                                 "magenta","orange","yellow","purple"),
                               mypch=15:19,mycex=2:10)
{

  Yg<-ergmm.fit$model$Yg
  Ym <- sociomatrix(Yg)
  n<-network.size(Yg)


  if(!missing(label) && length(label)==1 && is.character(label)){
    trycol <- unlist(get.vertex.attribute(Yg,label))
    if(!all(is.na(trycol))){
      label <- trycol
    }
  }
  
  if(!missing(vertex.pch) && length(vertex.pch)==1 && is.character(vertex.pch)){
    trycol <- as.numeric(as.factor(unlist(get.vertex.attribute(Yg,vertex.pch))))
    if(!all(is.na(trycol))){
      vertex.pch <- mypch[trycol]
    }
  }
  
  if(!missing(vertex.cex) && length(vertex.cex)==1 && is.character(vertex.cex)){
    trycol <- as.numeric(as.factor(unlist(get.vertex.attribute(Yg,vertex.cex))))
    if(!all(is.na(trycol))){
      vertex.cex <- mycex[trycol]
    }
  }


  if (!is.latent(ergmm.fit)) {
    ##
    ##   So regular (non-latent) call
    ##
    plot(ergmm.fit$samples$llk, type = "l",
         main = "MCMC log Likelihood Values", 
         xlab = "sample", ylab = "log likelihood",...)
  }
  else {
    ##
    ## So latent call
    ##
    if (missing(xlab)) 
      xlab <- ""
    if (missing(ylab)) 
      ylab <- ""
    if(mle){
      summ<-summary(ergmm.fit,point.est=c("mle","pmean"))
      z.pos <- summ$mle$Z
      if (missing(main)) 
        main <- paste("MLEs of Latent Positions of", 
                      deparse(substitute(ergmm.fit)))
    }
    else{
#### I have temporarily replaced MKL fits with posterior means. (PK)
##      z.pos <- ergmm.fit$Z.mkl
##      if (missing(main)) 
##        main <- paste("KL Latent Positions of", 
##                      deparse(substitute(ergmm.fit)))
      summ<-summary(ergmm.fit,point.est=c("pmean"))
      z.pos <- summ$pmean$Z
      if (missing(main)) 
        main <- paste("Posterior Mean Positions of", 
                      deparse(substitute(ergmm.fit)))
    }

    if(dim(z.pos)[2]==1){
      z.pos<-cbind(z.pos,z.pos)+matrix(rnorm(2*length(z.pos),0,sd(z.pos)/sqrt(length(z.pos)/2)),ncol=2)
    }

    if(missing(vertex.col)){
      if(is.latent.cluster(ergmm.fit))
        vertex.col <- cluster.col[summ$pmean$Z.K] 
    }
    else if(length(vertex.col)==1 && is.character(vertex.col)){
      trycol <- as.numeric(as.factor(unlist(get.vertex.attribute(Yg,vertex.col))))
      if(!all(is.na(trycol))){
        vertex.col <- cluster.col[trycol]
      }
    }
    
    if(density[1]){
      ##
      ##   Plot densities
      ##
      if(!require(KernSmooth,quietly=TRUE)){
       stop("The 'density' option requires the 'KernSmooth' package.")
      }
      if(density[1]>1)
      {
        ##
        ##   Plot 2-dimensional densities
        ##
        if(length(density)>1)
          opar <- par(mfrow=c(density[1],density[2]),mar=c(2.5,2.5,1,1))
        else
          opar <- par(mfrow=c(density,density),mar=c(2.5,2.5,1,1))
      }
      if(!only.subdens){
        Z.use <- cbind(as.vector(ergmm.fit$samples$Z[,,1]),as.vector(ergmm.fit$samples$Z[,,2]))
        plot(Z.use,type='n',xlab=xlab,ylab=ylab,xlim=range(Z.use),
             ylim=range(Z.use), ...)
        title(main=paste("Posterior density of",main), cex.main=0.7, ...)
        temp <- bkde2D(Z.use,0.2,c(201,201))
        image(temp$x1,temp$x2,temp$fhat,col=grey(seq(1,0,length=255)),add=TRUE)
        if(drawarrows)
          for (i in 1:n)
            for (j in 1:n)
              if (Ym[i,j])
                ergmm.midarrow(x0 = z.pos[i, 1], y0 = z.pos[i, 2], x1 = z.pos[j,1],
                               y1 = z.pos[j, 2], length = 0.1,col="yellow")
        box()
      }

      if(density[1]>1){
        Z.Z.K.v <- as.vector(t(ergmm.fit$Z.K))
        Z.proc.mean <- cbind(as.vector(ergmm.fit$samples$Z[,,1]),as.vector(ergmm.fit$samples$Z[,,2]))
        for(i in 1:ergmm.fit$ngroups){
          plot(Z.proc.mean,
               xlim=range(Z.proc.mean),ylim=range(Z.proc.mean),
               main=paste("Class",i),
               type='n', ...)
          temp <- bkde2D(Z.proc.mean[Z.Z.K.v==i,],0.2,c(101,101))
          use.col <- as.vector(col2rgb(cluster.col[i])/255)
          image(temp$x1,temp$x2,temp$fhat,add=TRUE,
                col=rgb(seq(1,use.col[1],length=255),
                  seq(1,use.col[2],length=255),
                  seq(1,use.col[3],length=255)))
          contour(temp$x1,temp$x2,temp$fhat,add=TRUE, nlevels=4,
                  drawlabels=FALSE,
                  col="white")
          box()
        }
        if(plotgraph)
          plot(ergmm.fit,label=label)
        par(opar)
      }
    }
    else if(contours>0){
      ## do a contours x contours array of contour plots of posterior densities
      opar <- par(mfrow=c(contours,contours),mar=c(0,0,1,0), omi=c(.5,.5,1,.5),
                  mgp=c(1.5,.25,0))
      temp.x.pos <- range(ergmm.fit$samples$Z[,,1])
      temp.y.pos <- range(ergmm.fit$samples$Z[,,2])
      if(!require(KernSmooth,quietly=TRUE)){
        stop("The 'contours' option requires the 'KernSmooth' package.")
      }
      for(k in 1:n){
        plot(x=ergmm.fit$samples$Z[,k,1],y=ergmm.fit$samples$Z[,k,2],
             xlab="",ylab="",type="n", asp=1,
             xlim=temp.x.pos,ylim=temp.y.pos,
             xaxt='n',yaxt='n')

        if((k-(contours^2)*trunc((k-1)/(contours^2)))/contours>(contours-1))
          axis(side=1)
        if(n-k < contours)
          axis(side=1)
        if(contours*trunc(k/contours)+1==k){
          axis(side=2)
        }
        est <- bkde2D(x=t(ergmm.fit$samples$Z[,k,]), bandwidth=c(0.2,0.2))
        est$fhat <- est$fhat/max(est$fhat,na.rm=TRUE)
        contour(est$x1, est$x2, est$fhat, add=TRUE, nlevels=5,
                drawlabels=FALSE,col=contour.color)
        for (i in 1:n)
          for (j in 1:n)
            if (Ym[i,j]){
              ergmm.midarrow(x0 = z.pos[i, 1], y0 = z.pos[i, 2], x1 = z.pos[j,1],
                             y1 = z.pos[j, 2], length = al/10,lwd=edge.lwd,col=edge.col)
            }
        if(!is.null(label) && length(label)==1 && is.character(label)){
          trycol <- unlist(get.vertex.attribute(ergmm.fit,label))
          if(!any(is.na(trycol))){
            label <- trycol
          }
        }
        if(is.null(label) | length(label)==1){
          label <- 1:n
        }
        text(z.pos[, 1], z.pos[, 2], label, col = label.col, cex = label.cex)
      }
      par(opar)
    }
    else{
      ##
      ##    If not densities or contours then
      ##    just plot the points with their MKL classes
      ##
      if (!is.null(comp.mat)) {
        procr <- function(Z, Zo) {
          A <- t(Z) %*% Zo %*% t(Zo) %*% Z
          A.eig <- eigen(A, symmetric = TRUE)
          A.sqrt <- A.eig$vec %*% diag(1/sqrt(A.eig$val)) %*% 
            t(A.eig$vec)
          Z %*% A.sqrt %*% t(Z) %*% Zo
        }
        z.pos <- procr(z.pos, comp.mat)
      }
      plot(z.pos, xlab = xlab, ylab = ylab,
           type = "n", asp = 1, main = main, ...)
      if(print.formula){
        aaa <- ergmm.fit$model$formula   
        xformula <- paste(aaa[2],aaa[1],aaa[-c(1:2)],collapse=" ")
        title(main = xformula, line = 1, cex.main = 0.7)
      }
      for (i in 1:n)
        for (j in 1:n)
          if (Ym[i,j]){
            ergmm.midarrow(x0 = z.pos[i, 1], y0 = z.pos[i, 2], x1 = z.pos[j,1],
                           y1 = z.pos[j, 2], length = al,lwd=edge.lwd,col=edge.col)
          }
      if (is.null(label)) 
        label <- 1:n
      if((!is.null(ergmm.fit$model$cluster)) & !pie){
        points(z.pos,cex=vertex.cex,col=vertex.col,pch=vertex.pch)
        if(label.col=="black") label.col <- "white"
      }
      else if(!pie){
        if(vertex.pch==19)
          points(z.pos,pch=19,cex=vertex.cex,col="white")
        points(z.pos,pch=vertex.pch,cex=vertex.cex,col=vertex.col)
        if(vertex.pch==19)
          points(z.pos,pch=21,cex=vertex.cex,col="black")
      }

      if(pie)
        for(i in 1:n){
          ergmm.drawpie(z.pos[i,],piesize,summ$pmean$Z.pZK[i,],n=50,cols=cluster.col)
          text(z.pos[, 1], z.pos[, 2], label, col = label.col, cex = label.cex)
        }
      else      
        text(z.pos[, 1], z.pos[, 2], label, col = label.col, cex = label.cex)
    }
  }
  invisible(NULL)
}
