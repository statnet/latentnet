library(latentnet)

n<-20

y.var<-1/16

y<-as.network(matrix(1,n,n),dir=TRUE)

cat("Euclidean:\n")

Z<-rbind(cbind(rnorm(n/2,3),rnorm(n/2,3)),
         cbind(rnorm(n/2,-3),rnorm(n/2,-3)))

dm<--as.matrix(dist(Z))

ym<-rnorm(n*n,2+dm,sqrt(y.var))

set.edge.value(y,"v",ym)
image(getYm(y,"v"))
y.fit<-ergmm(y~euclidean(d=2,G=2),response="v",family="normal",fam.par=list(var=y.var),
             verbose=TRUE)

Z.mkl<-plot(y.fit,Z.ref=Z)
points(Z,pch=5)
cat("Mean squared difference:",sum((Z.mkl-Z)^2),"\n")


cat("Bilinear:\n")

dm<-tcrossprod(Z)

ym<-rnorm(n*n,2+dm,sqrt(y.var))

set.edge.value(y,"v",ym)
image(getYm(y,"v"))
y.fit<-ergmm(y~bilinear(d=2,G=2),response="v",family="normal",fam.par=list(var=y.var),
             verbose=TRUE)

Z.mkl<-plot(y.fit,Z.ref=Z)
points(Z,pch=5)
cat("Mean squared difference:",sum((Z.mkl-Z)^2),"\n")
