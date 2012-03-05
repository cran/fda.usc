S.NW<-function (tt, h=NULL, Ker = Ker.norm,w=NULL,cv=FALSE) {
 if (is.matrix(tt)) {
    if (ncol(tt)!=nrow(tt)) {
      if (ncol(tt)==1) {
         tt=as.vector(tt)
         tt=abs(outer(tt,tt, "-"))}
      #else stop("Error: incorrect arguments passed")
    }}
 else if (is.vector(tt))    tt=abs(outer(tt,tt, "-"))
 else stop("Error: incorrect arguments passed")
      if (is.null(h)) {
    h=quantile(tt,probs=0.15,na.rm=TRUE)
    cat("h=");print(h)
    }
  if (cv)  diag(tt)=Inf
  k=Ker(tt/h)
  if (is.null(w)) w<-rep(1,nrow(k))
  k1<-sweep(k,1,w,FUN="*")   #antes un 2, aviso en prediccion
  S =k1/apply(k1,1,sum)
return(S)
}


