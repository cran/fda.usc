depth.mode=function(fdataobj,trim=0.25,metric=metric.lp,h=0.15,scale=TRUE,draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
nas<-apply(fdataobj$data,1,count.na)
if (any(nas))  {
   fdataobj$data<-fdataobj$data[!nas,]
   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
   }
data<-fdataobj[["data"]]
names1<-names2<-names<-fdataobj[["names"]]
names1$main<-"depth.mode median"
names2$main<-paste("depth.mode trim ",trim*100,"\u0025",sep="")
n<-nrow(data)
m<-ncol(data)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
tt=fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
mdist<-matrix(0.0,nrow=n,ncol=n)
mdist=metric(fdataobj,fdataobj,...)
  #   h<-max(mdist)*h
h=quantile(mdist+diag(Inf,nrow(mdist)),probs=h,na.rm=TRUE)
#    ans<-apply(mdist/h,1,skernel.norm) #see Kernel
ans<-Ker.norm(mdist/h)
ans<-apply(ans,1,sum,na.rm=TRUE)
if (scale) ans=as.vector(scale(ans,center=min(ans,na.rm=TRUE),scale=(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))))
k=which.max(ans)
med=data[k,]
lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
mtrim=apply(data[lista,],2,mean)
tr<-paste("mode.tr",trim*100,"\u0025",sep="")
med<-fdata(med,tt,rtt,names1)
mtrim<-fdata(mtrim,tt,rtt,names2)

rownames(med$data)<-"mode.med"
rownames(mtrim$data)<-tr
if (draw){
   ind1<-!is.nan(ans)
   ans[is.nan(ans)]=NA
   cgray=1-(ans-min(ans,na.rm=TRUE))/(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))
   plot(fdataobj[ind1,],col=gray(cgray[ind1]),main="mode Depth")
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"))
 }
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,"dep"=ans,"dist"=mdist)))
}

