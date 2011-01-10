
fregre.basis=function(fdataobj,y,basis.x=NULL,basis.b=NULL,lambda=0,
Lfdobj=vec2Lfd(c(0,0),rtt),...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
  call<-match.call()
  n = nrow(x)
  np <- ncol(x)
  if (n != (length(y))) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(rownames(x)))        rownames(x) <- 1:n
  if (is.null(colnames(x)))        colnames(x) <- 1:np
	if (is.null(basis.x))  {
	          nbasis.x=floor(np/6)
            basis.x=create.bspline.basis(rangeval=rtt,nbasis=nbasis.x,...)
            }
	if (is.null(basis.b))  {
     nbasis.b=floor(np/10)
     basis.b=create.bspline.basis(rangeval=rtt,nbasis=nbasis.b)
  }
	xmean=apply(x,2,mean)
  xcen=sweep(x,2,xmean,FUN="-")
	ymean=mean(y)
  ycen=y-ymean
  x.fd=Data2fd(argvals=tt,y=t(xcen),basisobj=basis.x)
  C=t(x.fd$coefs)
  J=inprod(basis.x,basis.b)
  vfunc="x."
  Z<-C%*%J
  Z=cbind(rep(1,len=n),Z)
  colnames(Z)<-1:ncol(Z)
  colnames(Z)[2:ncol(Z)]= paste(vfunc,substr(basis.b$type, 1,3), 1:basis.b$nbasis, sep = "")
if (lambda==0) {
       S=Z%*%solve(t(Z)%*%Z)%*%t(Z)
       yp2=S%*%y
       response="y"
       pf <- paste(response, "~", sep = "")
       for ( i in 1:length(colnames(Z))) pf <- paste(pf, "+", colnames(Z)[i], sep = "")
       object.lm=lm(formula=pf,data=data.frame(y,Z),x=TRUE,y=TRUE,...)
       yp=object.lm$fitted.values
       e<-object.lm$residuals
       b.est<-object.lm$coefficients[-1]
       beta.est=fd(b.est,basis.b)
       a.est<-object.lm$coefficients[1]
       df=basis.b$nbasis+1
	     sr2=sum(e^2)/(n-df)
	     r2=1-sum(e^2)/sum(ycen^2)
	     coefficients<-object.lm$coefficients
}
else {
       R=diag(0,ncol= basis.b$nbasis+1,nrow=basis.b$nbasis+1)
#       R[-1,-1]<-getbasispenalty(basis.b,Lfdobj) ############
       R[-1,-1]<-eval.penalty(basis.b,Lfdobj)
       Sb=t(Z)%*%Z+lambda*R
       eigchk(Sb)
       Cinv<-solve(Sb)
       Sb2=Cinv%*%t(Z)
       DD<-t(Z)%*%y
       S=Z%*%Sb2
       yp=S%*%y
       b.est=Sb2%*%y
       bet<-Cinv%*%DD
       a.est=b.est[1,1]
       #beta.est2=fd(b.est2[-1,1]*diff(rtt),basis.b)
       beta.est=fd(b.est[-1,1],basis.b)
       e=y-yp
       df=basis.b$nbasis+1
       sr2=sum(e^2)/(n-df)
       r2=1-sum(e^2)/sum(ycen^2)
       object.lm=list()
       object.lm$coefficients<-drop(b.est)
       object.lm$residuals<-drop(e)
       object.lm$fitted.values<-yp
       object.lm$y<-y
       object.lm$rank<-df
       object.lm$df.residual<-n-df
       colnames(Z)[1]="(Intercept)"
       vcov2=sr2*Cinv
       std.error=sqrt(diag(vcov2))
       t.value=b.est/std.error

       p.value= 2 * pt(abs(t.value),n-df, lower.tail = FALSE)
       coefficients<-cbind(b.est,std.error,t.value,p.value)
       colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
       class(object.lm)<-"lm"
b.est=b.est[-1]
names(b.est)<-rownames(coefficients)[-1]
     }
hat<-diag(hat(Z, intercept = TRUE),ncol=n)
out<-list("call"=call,"b.est"=b.est,"a.est"=a.est,"fitted.values"=yp,"H"=S,
  "residuals"=e,"df"=df,"r2"=r2,"sr2"=sr2,"y"=y,"fdataobj"=fdataobj,
  x.fd=x.fd,"beta.est"=beta.est,"basis.x.opt"=basis.x,"basis.b.opt"=basis.b,
  "J"=J,"lambda.opt"=lambda,lm=object.lm,coefficients=coefficients)
 class(out)="fregre.fd"
return(out)
}

