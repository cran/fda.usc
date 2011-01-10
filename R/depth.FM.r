depth.FM=function(fdataobj,trim=0.25,xeps=0.00000001,draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
data<-fdataobj[["data"]]
n<-nrow(data)
m<-ncol(data)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
t=fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
names1<-names2<-names<-fdataobj[["names"]]
names1$main<-"depht.FM median"
names2$main<-paste("FM trim ",trim*100,"\u0025",sep="")
d<-matrix(NA,nrow=n,ncol=m)
mdata<-data
for (i in 1:m){
       Fn=ecdf(mdata[,i])
#        d[,i]=Fn(fdata[,i])*(1-Fn(fdata[,i]-xeps))
         d[,i]=1-abs(0.5-Fn(mdata[,i]))
#        d[,i]<-1-abs(.5-rank(fdata[,i],ties.method="average")/n)
    }
    tr<-paste("FM.tr",trim*100,"\u0025",sep="")
    ans<-apply(d,1,mean)
#    rid<-rank(ans,ties.method="first")
    k=which.max(ans)
    med=data[k,]
    lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
    mtrim=apply(data[lista,],2,mean)
       if (draw){cgray=1-(ans-min(ans))/(max(ans)-min(ans))
       if (m==2){
   plot(range(data[,1]),range(data[,2]),type="n")
   text(data[,1],data[,2],round(ans,3),cex=0.75)
   points(rbind(mtrim),pch=19,col=gray(2*trim),cex=2)
   points(rbind(med),col=3,pch=20,cex=2)} else {
   plot(range(t),range(data),type="n",xlab="t",ylab="X(t)",main="FM Depth")
   for (i in 1:n) {lines(t,data[i,],col=gray(cgray[i]))}
   lines(t,mtrim,lwd=2,col="yellow")
   lines(t,med,col="red",lwd=2) }
   legend("topleft",legend=c(tr,"Median"),
    lwd=2,col=c("yellow","red"))
   }
med<-fdata(med,t,rtt,names1)
mtrim<-fdata(mtrim,t,rtt,names2)
rownames(med$data)<-"FM.med"
rownames(mtrim$data)<-tr
    return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,"dep"=ans)))
}


