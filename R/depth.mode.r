depth.mode=function(fdataobj,trim=0.25,metric=metric.lp,h=0.15,scale=TRUE,draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
data<-fdataobj[["data"]]
names1<-names2<-names<-fdataobj[["names"]]
names1$main<-"depth.mode median"
names2$main<-paste("depth.mode trim ",trim*100,"\u0025",sep="")
n<-nrow(data)
m<-ncol(data)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
t=fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
mdist<-matrix(0.0,nrow=n,ncol=n)
mdist=metric(fdataobj,...)
  #   h<-max(mdist)*h
h=quantile(mdist+diag(Inf,nrow(mdist)),probs=h,na.rm=TRUE)
#    ans<-apply(mdist/h,1,skernel.norm) #see Kernel
ans<-Ker.norm(mdist/h)
ans<-apply(ans,1,sum)
if (scale) ans=as.vector(scale(ans,center=min(ans),scale=(max(ans)-min(ans))))
k=which.max(ans)
med=data[k,]
lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
mtrim=apply(data[lista,],2,mean)
tr<-paste("mode.tr",trim*100,"\u0025",sep="")
if (draw){
       cgray=1-(ans-min(ans))/(max(ans)-min(ans))
    if (m==2){
   plot(range(data[,1]),range(data[,2]),type="n")
   text(data[,1],data[,2],round(ans,3),cex=0.75)
   points(rbind(mtrim),pch=19,col=gray(2*trim),cex=2)
   points(rbind(med),col=3,pch=20,cex=2)}
   else {
   plot(range(t),range(data),type="n",xlab="t",ylab="X(t)",main="Mode Depth")
   for (i in 1:n) {
   lines(t,data[i,],col=gray(cgray[i]))
   }
   lines(t,mtrim,lwd=2,col="yellow")
   lines(t,med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),
   lwd=2,col=c("yellow","red"))
   }
   }
med<-fdata(med,t,rtt,names1)
mtrim<-fdata(mtrim,t,rtt,names2)
rownames(med$data)<-"mode.med"
rownames(mtrim$data)<-tr
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,"dep"=ans,"dist"=mdist)))
    }




