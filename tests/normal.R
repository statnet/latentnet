library(latentnet)

n<-20

y<-as.network(matrix(1,n,n),dir=TRUE)

Z<-rbind(cbind(rnorm(n/2,4),rnorm(n/2,4)),
         cbind(rnorm(n/2,-4),rnorm(n/2,-4)))

plot(Z)

dm<-as.matrix(dist(Z))

ym<-rnorm(n*n,2-dm)

set.edge.value(y,"v",ym)
image(getYm(y,"v"))
y.fit<-ergmm(y~latent(d=2,G=2),response="v",family="normal",fam.par=list(var=1),
             verbose=TRUE)

plot(y.fit)
