S.KNN<-function (tt,knn=NULL,Ker=Ker.unif,cv=FALSE)      {
if (is.matrix(tt)) {
    if (ncol(tt)!=nrow(tt)) {
      if (ncol(tt)==1) {
         tt=as.vector(tt)
         tt=abs(outer(tt,tt, "-"))}
#      else stop("Error: incorrect arguments passed")
    }}
 else if (is.vector(tt))    tt=abs(outer(tt,tt, "-"))
 else stop("Error: incorrect arguments passed")
numgr=nrow(tt)
if (is.null(knn)) knn=floor(quantile(1:numgr,probs=0.05,na.rm=TRUE,type=4))
else if (knn<=0) stop("Error: incorrect knn value")
tol=1e-19
tol=diff(range(tt)*tol)
vec=apply(tt,1,quantile,probs=((knn+1)/numgr),type=4)+tol
rr=sweep(tt,1,vec,"/")
rr=Ker(rr)
if (cv) diag(rr)=0
res=rr/apply(rr,1,sum,na.rm=TRUE)
return(res)}

