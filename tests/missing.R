library(latentnet)
data(sampson)
ym<-as.matrix(samplike)
ym[1,5]<-NA
yg<-as.network(ym,matrix.type="adjacency")
samp.fit<-ergmm(yg~latent(d=2))

print(summary(samp.fit))
plot(samp.fit,labels=TRUE)
