depth.RPD=function(fdataobj,nproj=50,deriv=c(0,1),trim=0.25,dfunc2=depth.mode,
method="bspline",draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
nas<-apply(fdataobj$data,1,count.na)
if (any(nas))  {
   fdataobj$data<-fdataobj$data[!nas,]
   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
   }
data<-fdataobj[["data"]]
names1<-names2<-names<-fdataobj[["names"]]
names1$main<-"depth.RPD median"
names2$main<-paste("RPD trim ",trim*100,"\u0025",sep="")
n<-nrow(data)
m<-ncol(data)
modulo=function(z){sqrt(sum(z^2))}
if (is.null(n) || is.null(m)) stop("Input must be a matrix")
tt=fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
newfunc=array(NA,dim=c(n,m,length(deriv)))
for (ider in 1:length(deriv)){
  if (deriv[ider]==0) {newfunc[,,ider]=data}
  else {
  newfunc[,,ider]=fdata.deriv(data,nderiv=(ider-1),method=method,...)[["data"]]
  }
}
dep=rep(0.0,n)
vproject=matrix(0.0,nrow=n,ncol=length(deriv))
z=rnorm(m*nproj)
z=matrix(z,nrow=nproj,ncol=m)
modu=apply(z,1,modulo)
z=z/modu
pb=txtProgressBar(min=0,max=nproj,style=3)
for (j in 1:nproj){
    setTxtProgressBar(pb,j-0.5)
    for (ider in 1:length(deriv)){
        matriz=newfunc[,,ider]
        vproject[,ider]=matriz%*%z[j,]
         }

#    resul=dfunc2(vproject,...)
    par.dfunc=list()
    par.dfunc$fdataobj<-vproject
    par.dfunc$trim<-trim
    resul=do.call(dfunc2,par.dfunc)
    dep=dep+resul$dep
    setTxtProgressBar(pb,j)
        }
close(pb)
dep=dep/nproj
k=which.max(dep)
med=data[k,]
lista=which(dep>=quantile(dep,probs=trim))
mtrim=apply(data[lista,],2,mean,na.rm=TRUE)
tr<-paste("RPD.tr",trim*100,"\u0025",sep="")
med<-fdata(med,tt,rtt,names1)
mtrim<-fdata(mtrim,tt,rtt,names2)
rownames(med$data)<-"RPD.med"
rownames(mtrim$data)<-tr
if (draw){
   ans<-dep
   ind1<-!is.nan(ans)
   ans[is.nan(ans)]=NA
   cgray=1-(ans-min(ans,na.rm=TRUE))/(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))
   plot(fdataobj[ind1,],col=gray(cgray[ind1]),main="RPD Depth")
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"))
 }
return(invisible(list(median = med, lmed = k, mtrim = mtrim,
                      ltrim = lista, dep = dep, proj = z)))
}

