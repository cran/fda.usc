################################################################################
# dscale puede ser una funcion o un valor
################################################################################
depth.modep=function(lfdata,lfdataref=lfdata,h=NULL,metric,
                     par.metric=list(),method="euclidean",
                     scale=FALSE,trim=0.25,draw=FALSE,ask=FALSE) 
{
  equal<-identical(lfdata,lfdataref)
 if (class(lfdata)=="list"){
  lenl <- length(lfdata)
  lenl2 <- length(lfdataref)
  m0 <- nrow(lfdata[[1]])
  if (is.null(rownames(lfdata[[1]]$data))) 
    rownames(lfdata[[1]]$data) <- 1:m0
  nms <- rownames(lfdata[[1]]$data)
  nas <- NULL
  for (i in 1:lenl) {
    nas <- c(nas, na.action(na.omit(lfdata[[i]])))
  }
  nas <- unique(nas)
  nullans <- !is.null(nas)
  nam1<-names(lfdata)
  if (missing(metric)){
    if (is.fdata(lfdata[[nam1[1]]])) metric<-rep("metric.lp",len=lenl)
    else    metric<-rep("metric.dist",len=lenl)   
  }
  mdist2<-metric.ldata(lfdata=lfdata,lfdataref=lfdataref,metric=metric,par.metric=par.metric,method=method)
  if (equal) mdist<-mdist2
  else mdist<-metric.ldata(lfdata=lfdataref,lfdataref=lfdataref,metric=metric,par.metric=par.metric,method=method)
  m<-m0
  n<-nrow(lfdata[[1]])
}
else stop("Error in lfdata argument")
#attr(mdist, "par.metric") <-attr(atr, "par.metric")                                  
#attr(mdist, "call") <-attr(atr, "call")                                  
#attr(mdist, "method") <-method                                    
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
#if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
#h="h=0.15"
if (is.null(h))   {
  h<-0.1
  hq2=quantile(mdist,probs=h,na.rm=TRUE)
#print("es nulo h")  
}
else {
  #cat("no es nulo h ",h,"\n")    
  if (is.numeric(h))    hq2<-h  
  else hq2=quantile(mdist,probs=as.numeric(h),na.rm=TRUE)
}
#print(hq2)
class(mdist2)<-class(mdist)<-c("matrix","fdist")
dep<-Ker.norm(mdist2/hq2)    #### en (nxm) matrix (new X ref)
#print(ans)
dep<-apply(dep,1,sum,na.rm=TRUE)                                    
#print(ans)
#print(ans[1:3])

if (scale)   {
   dep2<-Ker.norm(mdist/hq2)
   dep2<-apply(dep2,1,sum,na.rm=TRUE)
   mn<-min(dep2,na.rm=TRUE)
   mx<-max(dep2,na.rm=TRUE)
#   scl<-mx-mn
#   ans=as.vector(scale(ans,center=mn,scale=scl))
#   ans2=as.vector(scale(ans2,center=mn,scale=scl))
   dep=as.vector(dep/mx)   
}
if (nullans) {
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-dep
        dep<-ans1      
        }
	names(dep)=nms
    k = which.max(dep) 
	nl=length(trim)
	tr<-paste("modep.tr",round(trim*100,2),"\u0025",sep="")
	lista=vector("list",nl)
	med=vector("list",lenl)
	mtrim=vector("list",lenl)
   for (ik in 1:lenl){
   med[[ik]]=lfdata[[ik]][k]
   mtrim[[ik]]=fdata(matrix(NA,ncol=length(lfdata[[ik]]$argvals),nrow=length(trim)),lfdata[[ik]]$argvals,lfdata[[ik]]$rangeval,lfdata[[ik]]$names)
   for (j in 1:nl){
   lista[[j]]=which(dep>=quantile(dep,probs=trim[j],na.rm=TRUE))
   	if (length(lista[[j]])==1) {
			mtrim[[ik]]$data[j,]<-lfdata[[ik]][lista[[j]]]$data
	if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
	}
	else   mtrim[[ik]]$data[j,]=apply(lfdata[[ik]]$data[lista[[j]],,drop=FALSE],2,mean,na.rm=TRUE)
	}
	rownames(med[[ik]]$data)="modep.med"
	rownames(mtrim[[ik]]$data)=tr
}
   if (nl>1) names(lista)=paste0("tr",trim*100)	
    if (draw) {
	  mf=5
		if (lenl>4) ask=TRUE
		if (ask) {par(mfrow = c(1, 1))
            dev.interactive()
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))}
		else{    mf<-switch(lenl,
			"1"={c(1,1)},
			"2"={c(1,2)},
			"3"={c(1,3)},
			"4"={c(2,2)})            
			par(mfrow =mf)                    
		}
        ans <- dep
        ind1 <- !is.nan(ans)
        ans[is.nan(ans)] = NA
		color=colorRampPalette(c("red","blue"))(nl+1)
        cgray = 1 - (ans - min(ans, na.rm = TRUE))/(max(ans,
            na.rm = TRUE) - min(ans, na.rm = TRUE))
for (ik in 1:lenl) {
   plot(lfdata[[ik]][ind1, ], col =gray(cgray[ind1]),lty=1, main = paste(nam1[ik],
  " modep Depth",sep=""))
   lines(mtrim[[ik]],lwd=3,col=color[-1],lty=1)
   lines(med[[ik]],col=color[1],lwd=3)
   legend("topleft",legend=c("Median",tr),lwd=3,col=color,box.col=0)
		}
    }

	out<- list(median = med, lmed = k, mtrim = mtrim,
	           ltrim = if (nl==1) unlist(lista) else lista, dep = dep,
	           "metric"=metric,"par.metric"=par.metric,"mdist"=mdist,"hq"=hq2)
	class(out)="depth"
    return(invisible(out))
} 


##################################################################################
depth.mode=function(fdataobj,fdataori=fdataobj,trim=0.25,metric=metric.lp,h=NULL,scale=FALSE,
draw=FALSE,...){    
#if (is.fdata(fdataobj)) {
# fdat<-TRUE
# nas<-is.na.fdata(fdataobj)
#if (any(nas))  {
#   fdataobj$data<-fdataobj$data[!nas,]
#   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
#   }
 if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
 nms<-rownames(fdataobj$data)
 m0<-nrow(fdataobj)
 fdataobj2=fdataobj
 fdataobj<-na.omit.fdata(fdataobj)
 fdataori<-na.omit.fdata(fdataori) 
 nas<-na.action(fdataobj)
 nullans<-!is.null(nas) 
# data<-fdataobj[["data"]]
# data2<-fdataori[["data"]]
 names1<-names2<-fdataobj[["names"]]
 names1$main<-"depth.mode median"
 names2$main<-paste("depth.mode trim ",trim*100,"\u0025",sep="")
 tt=fdataobj[["argvals"]]
 rtt<-fdataobj[["rangeval"]]
#print("is fdata") 
#}
#else { stop("no fdata class object")
#        data<-fdataobj
#        data2<-fdataori
#        fdat<-FALSE   
#  }
n<-nrow(fdataobj)
m<-ncol(fdataobj)
m2<-ncol(fdataori)
n2<-nrow(fdataori)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.matrix(metric)) mdist=metric                    
else  mdist=metric(fdataori,fdataori,...)
class(mdist)<-"matrix"

if (n==n2 & m==m2) {
  equal<-all(fdataobj$data==fdataori$data)
  if (equal) mdist2<-mdist
  else mdist2<-metric(fdataobj,fdataori,...)}
else  mdist2<-metric(fdataobj,fdataori,...)
if (is.null(h))   {
  h<-0.15
  hq2=quantile(mdist,probs=h,na.rm=TRUE)
  }
else {
  if (is.numeric(h))    hq2<-h  
  else hq2=quantile(mdist,probs=as.numeric(h),na.rm=TRUE)
}
class(mdist)<-class(mdist2)<-c("matrix","fdist")
dep<-Ker.norm(mdist2/hq2)    ####
dep<-apply(dep,1,sum,na.rm=TRUE)                                    
if (scale)   {
	dep2=Ker.norm(mdist/hq2)
	dep2=apply(dep2,1,sum,na.rm=TRUE)
   mn<-min(dep2,na.rm=TRUE)
   mx<-max(dep2,na.rm=TRUE)
   dep=as.vector(dep/mx)   
}
if (nullans) {
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-dep
        dep<-ans1      
        }
names(dep)<-nms		
k=which.max(dep)
med=fdataobj[k]
   nl=length(trim)
   lista=vector("list",nl)
   tr<-paste("mode.tr",round(trim*100,2),"\u0025",sep="")
   if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
   mtrim=fdata(matrix(NA,ncol=length(tt),nrow=length(trim)),tt,rtt,names2)
   for (j in 1:nl){
   lista[[j]]=which(dep>=quantile(dep,probs=trim[j],na.rm=TRUE))
   	if (length(lista[[j]])==1) {
			mtrim$data[j,]<-fdataobj2[lista[[j]]]$data
	if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
	}
	else   mtrim$data[j,]=apply(fdataobj2$data[lista[[j]],,drop=FALSE],2,mean,na.rm=TRUE)
	}
   rownames(med$data)<-"mode.med"
   rownames(mtrim$data)<-tr
   out<-list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=if (nl==1) unlist(lista) else lista,
"dep"=dep,"hq"=hq2,name="Mode") 
if (scale) out$dscale=mx
   out$trim <- trim
   out$name <- "mode"
   out$fdataobj=fdataobj
   out$fdataori=fdataori
  class(out)="depth"
  if (draw){
		plot.depth(out)
  }
  return(invisible(out))
}
################################################################################
mdepth.LD=function(x,xx=x,metric=metric.dist,
                   h=NULL,scale=FALSE,...){    
 if (is.vector(x)){
   if (all(xx==x)) x<-xx<-matrix(x,ncol=1)#stop("One of x or xx must be a matrix")
   else   {
 	m2=ncol(xx)
	if (length(x)!=m2) stop("Length of x does not match with dimension of xx")
	x=matrix(x,ncol=m2)
   }
 }    
		m2=ncol(xx)
		n2=nrow(xx)
		n<-nrow(x)
		m<-ncol(x)
if (is.null(rownames(x)))  rownames(x)<-1:nrow(x)
 nms<-rownames(x)
 x<-na.omit(x)  
 xx<-na.omit(xx)
 nas<-na.action(x)
 nullans<-!is.null(nas) 
        d <- ncol(x)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.matrix(metric)) {mdist=metric}
else {  mdist=metric(xx,xx,...)  }
class(mdist)<-"matrix"

if (n==n2 & m==m2) {
  equal<-all(x==xx)
  if (equal) mdist2<-mdist
  else mdist2<-metric(x,xx,...)}
else  mdist2<-metric(x,xx,...)
if (is.null(h))   {
  h<-0.15
  hq2=quantile(mdist,probs=h,na.rm=TRUE)
}
else {
  if (is.numeric(h))    hq2<-h  
  else hq2=quantile(mdist,probs=as.numeric(h),na.rm=TRUE)  
}
class(mdist)<-class(mdist2)<-c("matrix")
ans<-Ker.norm(mdist2/hq2)    ####
ans<-apply(ans,1,sum,na.rm=TRUE)                                    
mx = scale
if (scale)   {
#  mn<-min(ans,na.rm=TRUE)
  ans2=Ker.norm(mdist/hq2)
  ans2=apply(ans2,1,sum,na.rm=TRUE)
  mx<-max(ans2,na.rm=TRUE)
  ans=as.vector(ans/mx)   
}                                
if  (nullans){
        ans1<-rep(NA,len=n)
        ans1[-nas] <-ans 
        ans<-ans1      
        }
names(ans)<-nms   
out<-list("dep" = ans,"hq"=hq2,dscale=mx,
          x=x,xx=xx,name="LD")
class(out)<-"mdepth"
return(invisible(out))
}
##########################################################################
##########################################################################
