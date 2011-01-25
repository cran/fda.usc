fregre.np<-function(fdataobj,y,h=NULL,Ker=AKer.norm,metric=metric.lp,type.S=S.NW,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)

nas<-apply(fdataobj$data,1,count.na)
nas.g<-is.na(y)
if (is.null(names(y))) names(y)<-1:length(y)
if (any(nas) & !any(nas.g)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
  y<-y[bb]
   }
else {
if (!any(nas) & any(nas.g)) {
   cat("Warning: ",sum(nas.g)," values of group with NA are omited \n")
   bb<-!nas.g
   fdataobj$data<-fdataobj$data[bb,]
     y<-y[bb]
   }
else {
if (any(nas) & any(nas.g))  {
   bb<-!nas & !nas.g
   cat("Warning: ",sum(!bb)," curves  and values of group with NA are omited \n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
}}

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
mdist=metric(fdataobj,fdataobj,...)
ke<-deparse(substitute(Ker))
ty<-deparse(substitute(type.S))
if (is.null(h)) h=h.default(fdataobj,probs=c(0.05,0.05),len=1,metric = mdist,Ker =ke,
 type.S =ty,...)
    H =type.S(mdist,h,Ker,cv=FALSE)
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

