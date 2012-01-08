################################################################################
#  Depth measures for multivariate data:
#  depth.SD:  a simplicial depth (PD) measure.
#  depth.HD:  a half-space depth (HD) measure based on random projections.
#  depth.PD:  a projection depth (PD) measure based on random projections.
#  depth.MhD: a Mahalanobis depth (PD) measure.
################################################################################
data(tecator)
dat<-tecator$y[,2:3]
par(mfrow=c(2,2))
depth.SD(x=dat,draw=TRUE)
depth.HD(x=dat,draw=TRUE)
depth.PD(x=dat,draw=TRUE)
depth.MhD(x=dat,draw=TRUE)

dev.new()
par(mfrow=c(2,2))
dat<-tecator$y[,2:3]
nam<-names(dat)
lenseq<-71
rg1<-range(dat[,1])
seq1<-seq(rg1[1],rg1[2],len=lenseq)
rg2<-range(dat[,2])
seq2<-seq(rg2[1],rg2[2],len=lenseq)
newdat<-expand.grid(seq1,seq2)
dep1<-depth.SD(x=newdat,xx=dat)
dep2<-depth.HD(x=newdat,xx=dat)
dep3<-depth.PD(x=newdat,xx=dat)
dep4<-depth.MhD(x=newdat,xx=dat)
image(seq1,seq2,matrix(dep1$dep,lenseq),xlab=nam[1], ylab=nam[2],main="Simplicial depth")
points(dat[,1],dat[,2])
image(seq1,seq2,matrix(dep2$dep,lenseq),xlab=nam[1], ylab=nam[2],main="Half-space depth")
points(dat[,1],dat[,2])
image(seq1,seq2,matrix(dep3$dep,lenseq),xlab=nam[1], ylab=nam[2],main="Projection depth")
points(dat[,1],dat[,2])
image(seq1,seq2,matrix(dep4$dep,lenseq),xlab=nam[1], ylab=nam[2],main="Mahalanobis depth")
points(dat[,1],dat[,2])

# ycat<-ifelse(tecator$y[,"Fat"]<15,1,4)
# points(dat[,1],dat[,2],col=ycat,pch=ycat)

 