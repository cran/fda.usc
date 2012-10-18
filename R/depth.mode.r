#incorporar  Ker en la llamada
depth.mode=function(fdataobj,fdataori=fdataobj,trim=0.25,metric=metric.lp,h=NULL,scale=FALSE,
draw=FALSE,...){    
if (is.fdata(fdataobj)) {
 fdat<-TRUE
 nas<-apply(fdataobj$data,1,count.na)
#if (any(nas))  {
#   fdataobj$data<-fdataobj$data[!nas,]
#   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
#   }
 data<-fdataobj[["data"]]
 data2<-fdataori[["data"]]
 names1<-names2<-names<-fdataobj[["names"]]
 names1$main<-"depth.mode median"
 names2$main<-paste("depth.mode trim ",trim*100,"\u0025",sep="")
 tt=fdataobj[["argvals"]]
 rtt<-fdataobj[["rangeval"]]
}
else { data<-fdataobj;     fdat<-FALSE     }

n<-nrow(data)
m<-ncol(data)
m2<-ncol(data2)
n2<-nrow(data2)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")

if (is.matrix(metric)) {mdist=mdist2=metric}
else {
  mdist=metric(fdataobj,fdataobj,...)
  mdist2=metric(fdataobj,fdataori,...)
}
#aaa<-fdataobj==fdataori
#print(aaa)
#print(dim(mdist));print(dim(mdist2))
#   h<-max(mdist)*h

#hq=quantile(mdist+diag(Inf,nrow(mdist)),probs=h,na.rm=TRUE)
class(mdist2)<-class(mdist)<-"matrix"
#hq=quantile(mdist+diag(Inf,nrow(mdist)),probs=h)
if (is.null(h)) {
  h<-0.15
  hq2=quantile(mdist2,probs=h)
  }
else   hq2<-h
class(mdist2)<-class(mdist)<-c("matrix","fdist")
##    ans<-apply(mdist/h,1,skernel.norm) #see Kernel
ans<-Ker.norm(mdist2/hq2)    ####
ans<-apply(ans,1,sum,na.rm=TRUE)
          
#ans2<-Ker.norm(mdist/hq)    ####
#ans2<-apply(ans2,1,sum,na.rm=TRUE)
                                                                             
if (scale) ans=as.vector(scale(ans,center=min(ans,na.rm=TRUE),scale=(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))))
k=which.max(ans)
med=data[k,]
lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
mtrim=apply(data[lista,],2,mean)
tr<-paste("mode.tr",trim*100,"\u0025",sep="")
if (fdat) {
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
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),box.col=0)
 }
 }
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,
"dep"=ans,"dist"=mdist2,"hq"=hq2)))
}

