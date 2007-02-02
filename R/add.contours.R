ergmm.add.contours <- function(x,mycol=1:10,nlevels=4,...)
{
  Z.Z.K.v <- as.vector(t(x$Z.K))
  Z.proc.mean <- cbind(as.vector(x$Z[,1,]),as.vector(x$Z[,2,]))
  zlim <- -Inf
  for(i in 1:x$ngroups)
    {
      temp <- bkde2D(Z.proc.mean[Z.Z.K.v==i,],0.2,c(101,101))
      zlim <- max(zlim,temp$fhat)
      contour(temp$x1,temp$x2,temp$fhat,add=TRUE, nlevels=nlevels,
              drawlabels=FALSE, col=mycol[i],...)
    }
  return(invisible(NULL))
}
