depth.FM=function(fdataobj,trim=0.25,xeps=0.00000001,draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
nas<-apply(fdataobj$data,1,count.na)
if (any(nas))  {
   fdataobj$data<-fdataobj$data[!nas,]
   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
   }
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
ans<-apply(d,1,mean,na.rm=TRUE)
#    rid<-rank(ans,ties.method="first")
k=which.max(ans)
med=data[k,]
lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
mtrim=apply(data[lista,],2,mean,na.rm=TRUE)
med<-fdata(med,t,rtt,names1)
mtrim<-fdata(mtrim,t,rtt,names2)
rownames(med$data)<-"FM.med"
rownames(mtrim$data)<-tr
if (draw){
   ind1<-!is.nan(ans)
   ans[is.nan(ans)]=NA
   cgray=1-(ans-min(ans,na.rm=TRUE))/(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))
   plot(fdataobj[ind1,],col=gray(cgray[ind1]),main="FM Depth")
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"))
 }
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,"dep"=ans)))
}

