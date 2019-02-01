depth.RP<-function(fdataobj,fdataori=fdataobj,trim=0.25,nproj=50,proj="vexponential",
dfunc="TD1",par.dfunc=list(),scale=FALSE,draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
if (!is.fdata(fdataori)) fdataori=fdata(fdataori)

# nas<-is.na.fdata(fdataobj)
# if (any(nas))  {
#    fdataobj$data<-fdataobj$data[!nas,]
#    cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
#    }

 if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
 nms<-rownames(fdataobj$data)
 m0<-nrow(fdataobj)
 fdataobj2<-fdataobj
 fdataobj<-na.omit.fdata(fdataobj)
 fdataori<-na.omit.fdata(fdataori)
 nas<-na.action(fdataobj)
 nullans<-!is.null(nas)
#data<-fdataobj[["data"]]
#data2<-fdataori[["data"]]
n<-nrow(fdataobj)
m<-ncol(fdataobj)
m2<-ncol(fdataori)
n2<-nrow(fdataori)
if (is.null(n) && is.null(m))  stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
tt=fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
names1<-names2<-fdataobj[["names"]]
names1$main<-"depth.RP median"
names2$main<-paste("RP trim",trim*100,"\u0025",sep="")

 #### new
 if (is.fdata(proj)) {
  nproj<-nrow(proj)
  if (fdataobj$argvals!=proj$argvals || ncol(fdataobj)!=ncol(proj)) stop("Error in proj dimension")
  z<-proj
  }
else {
z<-rproc2fdata(nproj,tt,sigma=proj,norm=TRUE,...)
}
##
 dep=matrix(0.0,ncol=nproj,nrow=n)
 dep2=matrix(0.0,ncol=nproj,nrow=n2)
if (dfunc %in% c("Liu1","TD1","FM1")) {
Fn<-list()
 for (j in 1:nproj){
		valor=inprod.fdata(fdataobj,z[j])
		valor2=inprod.fdata(fdataori,z[j])	
        Fn[[j]]=ecdf(valor2)
         par.dfunc$x<-valor
         par.dfunc$Fn<-Fn[[j]]
         dep[,j]<-do.call(dfunc,par.dfunc)
		if (scale) {par.dfunc$x=valor2;dep2[,j]<-do.call(dfunc,par.dfunc)}
 }
  dep=apply(dep,1,mean)   #dep/nproj
  }
else if (dfunc %in% c("MhD1","LD1")) {
 for (j in 1:nproj){
 		valor=inprod.fdata(fdataobj,z[j])
		valor2=inprod.fdata(fdataori,z[j])
         par.dfunc$x<-drop(valor)
         par.dfunc$xx<-drop(valor2)
         dep[,j]<-do.call(dfunc,par.dfunc)
		if (scale) {par.dfunc$x=valor2;dep2[,j]<-do.call(dfunc,par.dfunc)}		 
#         dep<-dep+dp
 }
  dep=apply(dep,1,mean)   #dep/nproj
  }
  if (scale){    dep2=apply(dep2,1,mean);dep<-dep/max(dep2)    }
   names(dep)<-nms
   k=which.max(dep)
   med=fdataobj[k]
   if (nullans) {
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-dep
        dep<-ans1
        }
   nl=length(trim)
   lista=vector("list",nl)
   tr<-paste("RP.tr",round(trim*100,2),"\u0025",sep="")
   if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
   mtrim=fdata(matrix(NA,ncol=length(tt),nrow=length(trim)),tt,rtt,names1)
   for (j in 1:nl){
   lista[[j]]=which(dep>=quantile(dep,probs=trim[j],na.rm=TRUE))
   	if (length(lista[[j]])==1) {
			mtrim$data[j,]<-fdataobj2[lista[[j]]]$data
	if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
	}
	else   mtrim$data[j,]=apply(fdataobj2$data[lista[[j]],,drop=FALSE],2,mean,na.rm=TRUE)
	}
   rownames(med$data)<-"RP.med"
   rownames(mtrim$data)<-tr
   out=list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=if (nl==1) unlist(lista) else lista, "dep"=dep,"proj"=z,dfunc=dfunc,par.dfunc=par.dfunc)
   out$trim <- trim
   out$name <- "RP"
   out$fdataobj=fdataobj
   out$fdataori=fdataori
   class(out)="depth"
   if (draw){
   plot.depth(out)
			}
return(invisible(out))
}

