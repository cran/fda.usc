fdata.bootstrap <-function(fdataobj,statistic=func.mean,alpha=0.05,nb=200,
smo=0.0,draw=FALSE,draw.control=NULL,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
data<-fdataobj[["data"]]
estmues<-statistic(fdataobj,...)
nr<-nrow(fdataobj)
nc<-ncol(fdataobj)
tt =fdataobj[["argvals"]]
rtt=fdataobj[["rangeval"]]
names=fdataobj[["names"]]
distboot<-matrix(NA,nrow=nb)
estboot<-matrix(NA,nrow=nb,ncol=nc)
pb=txtProgressBar(min=0,max=nb,style=3)
for (i in 1:nb){
  setTxtProgressBar(pb,i-0.5)
  bmuestra<-fdataobj[sample(1:nr,size=nr,replace=TRUE),]
  if (smo>0) {
bmuestra[["data"]]<-bmuestra[["data"]]+mvrnorm(n=nr,rep(0,nc),var(data)*smo)}
stat<-statistic(bmuestra,...)
estboot[i,]<-stat[["data"]]
setTxtProgressBar(pb,i)}
close(pb)
center<-estmues
for (i in 1:nb){
  aux<-fdata(estboot[i,],tt,rtt)
  distboot[i]<-metric.lp(center,aux,...)
  }
dist<-max(distboot[rank(distboot)<=floor((1-alpha)*nb)])
if (draw){
if (is.null(draw.control)) draw.control=list("col"=c("grey","blue","cyan"),"lty"=c(2,1,1),"lwd"=c(1,2,1))
 plot(fdataobj,main="Bootstrap replies",lwd=draw.control$lwd[1],
 lty=draw.control$lty[1],col=draw.control$col[1])
 for(i in 1:nb){
if(distboot[i]<=dist) lines(tt,estboot[i,],lwd=draw.control$lwd[3],lty=draw.control$lty[3],col=draw.control$col[3])
#else lines(tt,estboot[i,],lwd=draw.control$lwd[4],lty=draw.control$lty[4],
#col=draw.control$col[4])
  }
lines(estmues,lwd=draw.control$lwd[2],lty=draw.control$lty[2],col=draw.control$col[2])
legend(x=min(tt),y=0.99*max(data),legend=c("fdata",stat$names$main,"IN"),
 lty=draw.control$lty,lwd=draw.control$lwd,col=draw.control$col,cex=0.9,box.col=0)
 }
return(list("statistic"=estmues,"dband"= dist,"rep.dist"=distboot,
    "resample"=fdata(estboot,tt,rtt,names),fdataobj=fdataobj))
}
#"center"=center,

