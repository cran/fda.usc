depth.RP=function(fdataobj,trim=0.25,nproj=50,xeps=0.0000001,draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
data<-fdataobj[["data"]]
n<-nrow(data)
m<-ncol(data)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
t=fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
names1<-names2<-names<-fdataobj[["names"]]
names1$main<-"depth.RP median"
names2$main<-paste("RP trim",trim*100,"\u0025",sep="")
 dep=rep(0.0,n)
 for (j in 1:nproj){
        z=rnorm(m)
        modulo=sum(z^2)
        z=z/sqrt(modulo)
        valor=data%*%z
        Fn=ecdf(valor)
        dep=dep+(Fn(valor)*(1-Fn(valor-xeps)))
 }
   dep=dep/nproj
   k=which.max(dep)
   med=data[k,]
   nl=length(trim)
   mtrim=matrix(NA,nrow=nl,ncol=m)
   for (j in 1:length(trim)) {
                    lista=which(dep>=quantile(dep,probs=trim[j],na.rm=TRUE))
                    mtrim[j,]=apply(data[lista,],2,mean)
                        }
   tr<-paste("RP.tr",trim*100,"\u0025",sep="")
   if (draw){cgray=1-(dep-min(dep))/(max(dep)-min(dep))
   if (m==2) {
   plot(range(data[,1]),range(data[,2]),type="n")
   text(data[,1],data[,2],round(dep,3),cex=0.75)
#   points(fdata[,1],fdata[,2],col=2,pch=19,cex=2)
   points(med[1],med[2],col=3,pch=20,cex=2)
   points(mtrim[,1],mtrim[,2],pch=19,col=gray(2*trim),cex=2)} else{
   plot(range(t),range(data),type="n",xlab="argvals",ylab="X(t)",main="RP Depth")
   for (i in 1:n) {lines(t,data[i,],col=gray(cgray[i]))}
   lines(t,mtrim,lwd=2,col="yellow")
   lines(t,med,col="red",lwd=2) }
   legend("topleft",legend=c(tr,"Median"),
    lwd=2,col=c("yellow","red"))
   }
med<-fdata(med,t,rtt,names1)
mtrim<-fdata(mtrim,t,rtt,names2)
rownames(med$data)<-"RP.med"
rownames(mtrim$data)<-tr
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,"dep"=dep,"proj"=z)))
}




