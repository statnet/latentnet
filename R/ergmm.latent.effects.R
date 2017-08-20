#  File R/ergmm.latent.effects.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
### effect-specific functions

## negative Euclidean distance
#' @importFrom stats dist
latent.effect.negative.Euclidean<-function(Z){
  -as.matrix(dist(Z))
}

dlpY.dZ.negative.Euclidean<-function(Z,dlpY.deta){
  n<-dim(Z)[1]
  d<-dim(Z)[2]
  Z.invdist<- as.matrix(dist(Z))
  Z.invdist[Z.invdist==0]<-Inf
  Z.invdist<-1/Z.invdist
  dlpY.dZ<-matrix(0,n,d)
  for(k in 1:d){
    Z.normdiff.k<-sapply(1:n,function(j)
                         sapply(1:n,function(i)
                                Z[i,k]-Z[j,k]))*Z.invdist
    dlpY.dZ[,k]<-dlpY.dZ[,k]+
      -sapply(1:n,function(i) crossprod(Z.normdiff.k[i,],dlpY.deta[i,]+dlpY.deta[,i]))
  }
  dlpY.dZ
}


## negative squared Euclidean distance
latent.effect.negative.Euclidean2<-function(Z){
  -as.matrix(dist(Z))^2
}

dlpY.dZ.negative.Euclidean2<-function(Z,dlpY.deta){
  n<-dim(Z)[1]
  d<-dim(Z)[2]
  dlpY.dZ<-matrix(0,n,d)
  for(k in 1:d){
    Z.norm.k<-sapply(1:n,function(j)
                     sapply(1:n,function(i)
                            Z[i,k]-Z[j,k]))*2
    dlpY.dZ[,k]<-dlpY.dZ[,k]+
      -sapply(1:n,function(i) crossprod(Z.norm.k[i,],dlpY.deta[i,]+dlpY.deta[,i]))
  }
  dlpY.dZ
}

## dot product (bilinear)
latent.effect.bilinear<-function(Z){
  tcrossprod(Z)
}

dlpY.dZ.bilinear<-function(Z,dlpY.deta){
  n<-dim(Z)[1]
  d<-dim(Z)[2]
  dlpY.dZ<-matrix(0,n,d)
  for(k in 1:d){
    Z.other<-sapply(1:n,function(j)
                         sapply(1:n,function(i)
                                Z[j,k]))
    dlpY.dZ[,k]<-dlpY.dZ[,k]+
      sapply(1:n,function(i) crossprod(Z.other[i,],dlpY.deta[i,]+dlpY.deta[,i]))
  }
  dlpY.dZ
}

## Dispatcher and information functions

latent.effect.IDs<-list(negative.Euclidean=1,
                        bilinear=2,
                        negative.Euclidean2=3)

latent.effect.names<-c("negative.Euclidean",
                       "bilinear",
                       "negative.Euclidean2")

latent.effect.invariances <- list(negative.Euclidean = c("translation", "rotation", "reflection"),
                                  bilinear = c("rotation", "reflection"),
                                  negative.Euclidean = c("translation", "rotation", "reflection"))


latent.effect.fs<-c(latent.effect.negative.Euclidean,
                    latent.effect.bilinear,
                    latent.effect.negative.Euclidean2)

dlpY.dZ.fs<-c(dlpY.dZ.negative.Euclidean,
              dlpY.dZ.bilinear,
              dlpY.dZ.negative.Euclidean2)
