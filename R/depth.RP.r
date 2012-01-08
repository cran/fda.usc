depth.RP<-function(fdataobj,fdataori=fdataobj,trim=0.25,nproj=50,proj=1,xeps=0.0000001,draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
if (!is.fdata(fdataori)) fdataobj=fdata(fdataori) 
nas<-apply(fdataobj$data,1,count.na)
if (any(nas))  {
   fdataobj$data<-fdataobj$data[!nas,]
   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
   }
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
names1$main<-"depth.RP median"
names2$main<-paste("RP trim",trim*100,"\u0025",sep="")
 dep=rep(0.0,n)
 dep2=rep(0.0,n2)
 #### new
 if (is.fdata(proj)) {
  nproj<-nrow(proj)
  if (fdataobj$argvals!=proj$argvals || ncol(fdataobj)!=ncol(proj)) stop("Error en proj dimension")
  z<-proj
  }
else {
z<-rproc2fdata(nproj,t,sigma=proj,norm=TRUE,...)
}
##
# z<-matrix(NA,nproj,m)
Fn<-list()   
 for (j in 1:nproj){
    #    z[j,]=rnorm(m)
     #   modulo=sum(z[j]^2)
     #   z[j,]=z[j,]/sqrt(modulo)
        valor=data%*%z$data[j,]   
        valor2=data2%*%z$data[j,]   
        Fn[[j]]=ecdf(valor2)
        dep=dep+(Fn[[j]](valor)*(1-Fn[[j]](valor-xeps)))
        dep2=dep2+(Fn[[j]](valor2)*(1-Fn[[j]](valor2-xeps)))        
 }
   dep=dep/nproj
   dep2=dep2/nproj   
   k=which.max(dep2)
   med=data2[k,]
   nl=length(trim)
#   mtrim=matrix(NA,nrow=nl,ncol=m)
#   for (j in 1:length(trim)) {
#                    lista=which(dep2>=quantile(dep2,probs=trim[j],na.rm=TRUE))
#                    if (length(lista)==1) {
#                      mtrim[j,]<-data2[lista,]
#                      if (draw) {draw=FALSE;warning("The plot is not shown")}
#                    }
#                    else mtrim[j,]=apply(data2[lista,],2,mean,na.rm=TRUE)
 #                       }
  lista=which(dep2>=quantile(dep2,probs=trim,na.rm=TRUE))
   mtrim=apply(data2[lista,],2,mean,na.rm=TRUE)
   tr<-paste("RP.tr",trim*100,"\u0025",sep="")
   med<-fdata(med,t,rtt,names1)
   mtrim<-fdata(mtrim,t,rtt,names2)
   rownames(med$data)<-"RP.med"
   rownames(mtrim$data)<-tr
if (draw){
   ans<-dep2
   ind1<-!is.na(ans)
   cgray=1-(ans-min(ans,na.rm=TRUE))/(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))
   plot(fdataori[ind1,],col=gray(cgray[ind1]),main="RP Depth")
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),box.col=0)
 }
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,
"dep"=dep,"dep.ori"=dep2,"proj"=z,"Fn"=Fn)))
}
