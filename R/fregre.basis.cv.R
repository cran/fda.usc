fregre.basis.cv=function(fdataobj,y,basis.x=NULL,basis.b=NULL,type.basis=NULL,lambda=0,
Lfdobj=vec2Lfd(c(0,0),rtt),type.CV=GCV.S,par.CV=list(trim=0),...){
call<-match.call()
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
n = nrow(x)
np <- ncol(x)
if (n != (length(y))) stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(rownames(x)))        rownames(x) <- 1:n
if (is.null(colnames(x)))        colnames(x) <- 1:np
if (is.matrix(y)) y=as.vector(y)
if (is.null(basis.x))  {
  nbasis1=floor(np/6)
  if (!is.null(type.basis))  {
    aa1 <- paste("create.", type.basis[1], ".basis", sep = "")
    as <- list()
    as[[1]] <- range(tt)
    names(as)[[1]] <- "rangeval"
    as[[2]] <- nbasis1
    names(as)[[2]] <- "nbasis"
    basis.x=do.call(aa1, as)
    }
  else   basis.x=create.bspline.basis(rangeval=rtt,nbasis=nbasis1,...)
}
if (is.null(basis.b))  {
  nbasis2=floor(np/10)
	if (!is.null(type.basis))  {
    if (length(type.basis)>1)    aa1 <- paste("create.",type.basis[2], ".basis", sep = "")
    else     aa1 <- paste("create.",type.basis, ".basis", sep = "")
    as <- list()
    as[[1]] <- rtt
    names(as)[[1]] <- "rangeval"
    as[[2]] <- nbasis2
    names(as)[[2]] <- "nbasis"
    basis.b=do.call(aa1, as)
    }
  else basis.b=create.bspline.basis(rangeval=range(tt),nbasis=nbasis2)
}
lenlambda=length(lambda)
a1=list()
a2=list()
if (is.vector(basis.x))  {
 lenbasis.x=length(basis.x)
 for (nb.x in 1:lenbasis.x) {
   if (!is.null(type.basis))  {
     aa1 <- paste("create.", type.basis[1], ".basis", sep = "")
     as <- list()
     as[[1]] <- rtt
     names(as)[[1]] <- "rangeval"
     as[[2]] <- basis.x[nb.x]
     names(as)[[2]] <- "nbasis"
     a1[[nb.x]]=do.call(aa1, as)
  }
 else  a1[[nb.x]]=create.bspline.basis(rangeval=rtt,nbasis=basis.x[nb.x])
 }
 basis.x=a1
}
else {
  lenbasis.x=1
  basis.x=list(basis.x)
  }
if (is.vector(basis.b))  {
  lenbasis.y=length(basis.b)
  for (nb.y in 1:lenbasis.y) {
    if (!is.null(type.basis))  {
      if (length(type.basis)>1)       aa1 <- paste("create.", type.basis[2], ".basis", sep = "")
      else     aa1 <- paste("create.", type.basis[1], ".basis", sep = "")
      as <- list()
      as[[1]] <- rtt
      names(as)[[1]] <- "rangeval"
      as[[2]] <- basis.b[nb.y]
      names(as)[[2]] <- "nbasis"
      a2[[nb.y]]=do.call(aa1, as)
      }
  else  a2[[nb.y]]=create.bspline.basis(rangeval=rtt,nbasis=basis.b[nb.y],...)
 }
 basis.b=a2
}
else {
 lenbasis.y=1
 basis.b=list(basis.b)
}
gcv=array(Inf,dim=c(lenbasis.x,lenbasis.y,lenlambda))
pr=Inf
i.lambda.opt=1;i.nb.y.opt=1;i.nb.x.opt=1
 xx<-fdata.cen(fdataobj)
	xmean=xx[[2]]
  xcen=xx[[1]]
    ymean=mean(y)
  ycen=y-ymean
for (nb.x in 1:lenbasis.x) {
	x.fd=Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x[[nb.x]])
  C=t(x.fd$coefs)
  Cm=t(mean.fd(x.fd)$coefs)
  for (nb.y in 1:lenbasis.y) {
   	 J=inprod(basis.x[[nb.x]],basis.b[[nb.y]])
	   Z=C%*%J
     Z=cbind(rep(1,len=n),Z)
	   if (min(lambda)!=0) {
       R=diag(0,ncol= basis.b[[nb.y]]$nbasis+1,nrow=basis.b[[nb.y]]$nbasis+1)
       R[-1,-1]<-eval.penalty(basis.b[[nb.y]],Lfdobj)
                    }
     else R=0
     for (k in 1:lenlambda) {
       Sb=t(Z)%*%Z+lambda[k]*R
       eigchk(Sb)
       Cinv<-solve(Sb)
       Sb2=Cinv%*%t(Z)
       par.CV$S<-Z%*%Sb2
       par.CV$y<-y
       gcv[nb.x,nb.y,k]<- do.call(type.CV,par.CV)
     if (gcv[nb.x,nb.y,k]<pr) {
          pr=gcv[nb.x,nb.y,k]
          lambda.opt=lambda[k]
          basis.b.opt=basis.b[[nb.y]]
          basis.x.opt=basis.x[[nb.x]]
          Sb.opt=Sb2
          Z.opt=Z
          Cm.opt=Cm
          J.opt=J
          Cinv.opt=Cinv
     }    }
    }  }

    l = which.min(gcv)
    gcv.opt=min(gcv)
    S=Z.opt%*%Sb.opt
    DD<-t(Z.opt)%*%y
    yp=S%*%y
    b.est=Sb.opt%*%y
    bet<-Cinv.opt%*%DD
    rownames(b.est)<-1:nrow(b.est)
    rownames(b.est)[1]<- "(Intercept)"
    #beta.est2=fd(b.est2[-1,1]*diff(rtt),basis.b)
    beta.est=fd(b.est[-1,1],basis.b.opt)
    a.est=b.est[1,1]
    e=y-yp
    df=basis.b.opt$nbasis+1
    sr2=sum(e^2)/(n-df)
    r2=1-sum(e^2)/sum(ycen^2)
    object.lm=list()
    object.lm$coefficients<-drop(b.est)
    object.lm$residuals<-drop(e)
    object.lm$fitted.values<-yp
    object.lm$y<-y
    object.lm$rank<-df
    object.lm$df.residual<-n-df
    vfunc=call[[2]]
    colnames(Z.opt)<-1:ncol(Z.opt)
    colnames(Z.opt)[2:ncol(Z.opt)]= paste(vfunc,".",basis.b.opt$names, sep = "")
    colnames(Z.opt)[1]="(Intercept)"
    vcov2=sr2*Cinv.opt
    std.error=sqrt(diag(vcov2))
    t.value=b.est/std.error
    p.value= 2 * pt(abs(t.value),n-df, lower.tail = FALSE)
    coefficients<-cbind(b.est,std.error,t.value,p.value)
    colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    rownames(coefficients)<-colnames(Z.opt)
    class(object.lm)<-"lm"
b.est=b.est[-1]
names(b.est)<-rownames(coefficients)[-1]
####no es la matriz de hat ya que yp=a.est+S*y y no yp=S*y
#yp=mean(y)+S%*%ycen
#b.est=Sb.opt%*%y #idem beta.est$coefs
#beta.est=fd(b.est*diff(rtt),basis.b.opt)
#  a.est=ymean-Cm%*%J%*%beta.est$coefs
#a.est=ymean-Cm.opt%*%J.opt%*%b.est
#e=y-yp
#df=sum(diag(S))+1
#r=x.fd[[2]][[3]]
#xh=cbind(rep(1,len=nrow(Z.opt)),Z.opt)
#betah=c(a.est,b.est)
#colnames(xh)[1]="(Intercept)"
#vcov2=sr2*solve(t(xh)%*%xh)
#std.error=sqrt(diag(vcov2))
#t.value=betah/std.error
#p.value= 2 * pt(abs(t.value),n-df, lower.tail = FALSE)
#result<-cbind(betah,std.error,t.value,p.value)
#colnames(result) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
out<-list("call"=call,"b.est"=b.est,"a.est"=a.est,"fitted.values"=yp,"H"=S,
"residuals"=e,"df"=df,"r2"=r2,"sr2"=sr2,"y"=y,"fdataobj"=fdataobj,"gcv"=gcv,
"lambda.opt"=lambda.opt,"gcv.opt"=gcv.opt,"coefficients"=coefficients,
"basis.x.opt"=basis.x.opt,"basis.b.opt"=basis.b.opt,"J"=J.opt,"beta.est"=beta.est,
"lm"=object.lm,"mean"=xmean)
class(out)="fregre.fd"
return(out)
}

