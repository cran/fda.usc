depth.RPD=function(fdataobj,nproj=50,deriv=c(0,1),trim=0.25,dfunc2=depth.mode,
method="bspline",draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
data<-fdataobj[["data"]]
names1<-names2<-names<-fdataobj[["names"]]
names1$main<-"depth.RPD median"
names2$main<-paste("RPD trim ",trim*100,"\u0025",sep="")
n<-nrow(data)
m<-ncol(data)
modulo=function(z){sqrt(sum(z^2))}
if (is.null(n) || is.null(m)) stop("Input must be a matrix")
t=fdataobj[["argvals"]]
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
mtrim=apply(data[lista,],2,mean)
tr<-paste("RPD.tr",trim*100,"\u0025",sep="")
if (draw){
 dev.new()
 cgray=1-(dep-min(dep))/(max(dep)-min(dep))
 if (m==2){
  plot(range(data[,1]),range(data[,2]),type="n",
     xlab=colnames(data)[1],ylab=colnames(data)[2])
  points(data[,1],data[,2],col=gray(cgray))
  points(rbind(mtrim),pch=19,col=gray(2*trim),cex=2)
  points(rbind(med),col=3,pch=20,cex=2)}
else {
plot(range(t),range(data),type="n",xlab="t",ylab="X(t)",main="RPD Depth")
for (i in 1:n) {lines(t,data[i,],col=gray(cgray[i]))}
lines(t,mtrim,lwd=2,col="yellow")
lines(t,med,col="red",lwd=2)
legend("topleft",legend=c(tr,"Median"),
     lwd=2,col=c("yellow","red"))
 }}
med<-fdata(med,t,rtt,names1)
mtrim<-fdata(mtrim,t,rtt,names2)
rownames(med$data)<-"RPD.med"
rownames(mtrim$data)<-tr
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,"dep"=dep,"proj"=z)))
}
