################################################################################
################################################################################
depth.FM=function(fdataobj,fdataori=fdataobj,trim=0.25,xeps=0.00000001,draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
if (!is.fdata(fdataori)) fdataobj=fdata(fdataori)
#nas<-apply(fd-ataobj$data,1,count.na)
#if (any(nas))  {
#   fdataobj$data<-fdataobj$data[!nas,]
#   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
#   }                 
             
data<-fdataobj[["data"]]
data2<-fdataori[["data"]]
n<-nrow(data)
m<-ncol(data)
m2<-ncol(data2)
n2<-nrow(data2)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
t=fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
names1<-names2<-names<-fdataobj[["names"]]
names1$main<-"depht.FM median"
names2$main<-paste("FM trim ",trim*100,"\u0025",sep="")
d<-matrix(NA,nrow=n,ncol=m)
d2<-matrix(NA,nrow=n2,ncol=m)
mdata<-data
Fn<-list()
for (i in 1:m){
       Fn[[i]]=ecdf(data2[,i])
#        d[,i]=Fn(fdata[,i])*(1-Fn(fdata[,i]-xeps))
#        d[,i]<-1-abs(.5-rank(fdata[,i],ties.method="average")/n)
         d[,i]=1-abs(0.5-Fn[[i]](mdata[,i]))
         d2[,i]=1-abs(0.5-Fn[[i]](data2[,i]))  
    }
tr<-paste("FM.tr",trim*100,"\u0025",sep="")
ans<-apply(d,1,mean,na.rm=TRUE)[1:n]
ans2<-apply(d2,1,mean,na.rm=TRUE)
#    rid<-rank(ans,ties.method="first")
##################   
k=which.max(ans2)    
med=data2[k,]
lista=which(ans2>=quantile(ans2,probs=trim,na.rm=TRUE))
#if (length(lista)==1) {
#  mtrim<-data2[lista,]
#  if (draw) {draw=FALSE;warning("The plot is not shown")}
#  }
#else mtrim=apply(data2[lista,],2,mean,na.rm=TRUE)
mtrim=apply(data2[lista,],2,mean,na.rm=TRUE)
med<-fdata(med,t,rtt,names1)
mtrim<-fdata(mtrim,t,rtt,names2)
rownames(med$data)<-"FM.med"
rownames(mtrim$data)<-tr
if (draw){
   ind1<-!is.nan(ans2)
  ans2[is.nan(ans2)]=NA
   cgray=1-(ans2-min(ans2,na.rm=TRUE))/(max(ans2,na.rm=TRUE)-min(ans2,na.rm=TRUE))
#   plot(fdataori[ind1,],col=gray(cgray[ind1]),main="FM Depth")
   plot(fdataori,col=gray(cgray),main="FM Depth")
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),
   box.col=0)
 }
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,
"dep"=ans,"dep.ori"=ans2,"Fn"=Fn)))
} 

