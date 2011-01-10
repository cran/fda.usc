fregre.np=function(fdataobj,y,h=NULL,Ker=AKer.norm,metric=metric.lp,type.S=S.NW,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
C<-match.call()
mf <- match.call(expand.dots = FALSE)
m<-match(c("fdataobj", "y","h","Ker","metric","type.S"),names(mf),0L)
#    if (is.vector(x))         x <- t(as.matrix(x))
n = nrow(x)
np <- ncol(x)
if (n != (length(y)))    stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(rownames(x)))         rownames(x) <- 1:n
if (is.null(colnames(x)))         colnames(x) <- 1:np
mdist=metric(x,x,...)
if (m[6]==0){
   		if (is.null(h)) h=quantile(mdist,probs=0.05,na.rm=TRUE)
   		 H =type.S(mdist,h,Ker,cv=FALSE)
   		}
    else {
       if (C[[m[6]]]=="S.KNN"){
            if (is.null(h))  {
                h=floor(quantile(1:n,probs=0.05,na.rm=TRUE,type=4))
                H =type.S(mdist,h,Ker,cv=FALSE)
                }
            else {H =type.S(mdist,h,Ker,cv=FALSE)}
            }
       else {	if (is.null(h)) h=quantile(mdist,probs=0.05,na.rm=TRUE)
        H =type.S(mdist,h,Ker,cv=FALSE)}
            }
  	yp=H%*%y
    ycen=y-mean(y)
 	  e=matrix(y,ncol=1)-yp
    df=traza(H)
   	sr2=sum(e^2)/(n-df)
  	r2=1-sum(e^2)/sum(ycen^2)
out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df"=df,"r2"=r2,"sr2"=sr2,
"y"=y,"fdataobj"=fdataobj,"mdist"=mdist,"Ker"=Ker,"metric"=metric,"type.S"=type.S,"h.opt"=h,"m"=m)
class(out)="fregre.fd"
return(out)
}



