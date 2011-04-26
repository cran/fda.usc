pls.fdata<- function(fdataobj,y,ncomp=2,norm=TRUE,...){
if (!is.fdata(fdataobj)) stop("No fdata class")
 nas1<-apply(fdataobj$data,1,count.na)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
X<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
nam<-fdataobj[["names"]]
mm<-fdata.cen(fdataobj)
xmean<-mm$meanX
Xcen.fdata <- mm$Xcen
n <- nrow(Xcen.fdata);J <- ncol(Xcen.fdata)
Jmin<-min(c(ncomp,J,n))####
eigenres <- svd(Xcen.fdata$data)
    cnames<-paste("PLS",colnames(fdata.cen(fdataobj)[[1]]$data),sep="")
    ycen = y - mean(y)
    response = "y"
    df<-data.frame(y,Xcen.fdata$data)
    colnames(df)<-c("y",cnames)
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf,"+",cnames[i],sep="")
    res<-plsr(as.formula(pf),ncomp=ncomp,data=df,...)
class(res$loadings)<-"matrix"
vs<-fdata(t(res$loadings),tt,rtt,list(main="pls.fdata",xlab="t",ylab="rotationloadings"))
scores <- matrix(0,ncol=J,nrow=n)
if (norm) {
          no<-norm.fdata(vs)
          vs$data<-sweep(vs$data,1,drop(no),"/")
          scores[,1:Jmin] <-inprod.fdata(Xcen.fdata,vs)
          }
else     {scores[,1:Jmin] <-inprod.fdata(Xcen.fdata,vs)}
#else     {scores[,1:Jmin] <-Xcen.fdata$data%*%t(eigenres$v)*diff(rtt)/J}
colnames(scores)<-paste("PLS",1:J, sep = "")
return(list("rotation" = vs, "x" = scores,"res.pls"=res,"fdataobj.cen"=Xcen.fdata,"mean"=xmean))
}



