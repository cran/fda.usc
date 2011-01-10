pc.svd.fdata<- function(fdataobj,norm=TRUE){
if (!is.fdata(fdataobj)) stop("No fdata class")
X<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
nam<-fdataobj[["names"]]
Xcen <- sweep(X,2,apply(X,2,mean),"-")
Xcen.fdata<-fdata(Xcen,tt,rtt,nam)
n <- nrow(Xcen);J <- ncol(Xcen)
Jmin<-min(c(J,n))
eigenres <- svd(Xcen)
lambdas <- eigenres$d^2/(n-1)
vs<-fdata(t(eigenres$v),tt,rtt,list(main="pc.svd.fdata",xlab="t",ylab="rotation"))
scores <- matrix(0,ncol=J,nrow=n)
if (norm) {
          no<-norm.fdata(vs)     #así norma 1, con mean(norm.fdata(vs)) aprox 1
          vs$data<-sweep(vs$data,1,drop(no),"/")
          scores[,1:Jmin] <-inprod.fdata(Xcen.fdata,vs)
          }
else     {scores[,1:Jmin] <-inprod.fdata(Xcen.fdata,vs)}
colnames(scores)<-paste("PC",1:J, sep = "")
return(list("lambdas"=lambdas,"rotation"=vs,"x"=scores,"fdataobj.cen"=Xcen.fdata))
}
