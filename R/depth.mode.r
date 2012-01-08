#incorporar  Ker en la llamada
depth.mode=function(fdataobj,fdataori=fdataobj,trim=0.25,metric=metric.lp,h=NULL,scale=FALSE,
draw=FALSE,...){    
if (is.fdata(fdataobj)) {
 fdat<-TRUE
# nas<-apply(fdataobj$data,1,count.na)
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
  mdist=metric(fdataori,fdataori,...)
  mdist2=metric(fdataobj,fdataori,...)
}
aaa<-fdataobj==fdataori
#print(dim(mdist));print(dim(mdist2))
#   h<-max(mdist)*h   
#hq=quantile(mdist+diag(Inf,nrow(mdist)),probs=h,na.rm=TRUE)
class(mdist2)<-class(mdist)<-"matrix"
#hq=quantile(mdist+diag(Inf,nrow(mdist)),probs=h)
if (is.null(h)) {
  h<-0.15
  hq2=quantile(mdist,probs=h)
  }
else   hq2<-h
class(mdist2)<-class(mdist)<-c("matrix","fdist")
##    ans<-apply(mdist/h,1,skernel.norm) #see Kernel      
ans<-Ker.norm(mdist/hq2)    ####
ans<-apply(ans,1,sum,na.rm=TRUE)                                    
ans2<-Ker.norm(mdist2/hq2)    ####
ans2<-apply(ans2,1,sum,na.rm=TRUE)                                                                           
if (scale) 
  {
   mn<-min(ans,na.rm=TRUE)
   mx<-max(ans,na.rm=TRUE)
   scl<-mx-mn
   ans=as.vector(scale(ans,center=mn,scale=scl))
   ans2=as.vector(scale(ans2,center=mn,scale=scl))
   }
k=which.max(ans)
med=data2[k,]
lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
if (length(lista)==1) {
  mtrim<-data2[lista,]
  if (draw) {draw=FALSE;warning("The plot is not shown")}
  }
else mtrim=apply(data2[lista,],2,mean,na.rm=TRUE)
tr<-paste("mode.tr",trim*100,"\u0025",sep="")
if (fdat) {
med<-fdata(med,tt,rtt,names1)
mtrim<-fdata(mtrim,tt,rtt,names2)
rownames(med$data)<-"mode.med"
rownames(mtrim$data)<-tr
if (draw){
    if (!scale){
     mn<-min(ans,na.rm=TRUE)
     mx<-max(ans,na.rm=TRUE)
    scl<-mx-mn
    }     
   ind1<-!is.nan(ans)
   ans[is.nan(ans)]=NA
   cgray=1-(ans-mn)/(scl)
   plot(fdataori[ind1,],col=gray(cgray[ind1]),main="mode Depth")
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),box.col=0)
 }
 }
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,
"dep"=ans2,"dep.ori"=ans,"dist"=mdist2,"hq"=hq2)))
}
