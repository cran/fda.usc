pc.cor=function(out,y=NULL,l=NULL,draw=TRUE,...) {
 if (class(out)=="fregre.fd") {
    if (!is.null(l)) l=l
    else l=out$l
   y=out$y
   pr.com=out$svd.fdata
   }
 else {
   if (!is.fdata(out)) stop("No fdata class")
   if (is.null(l)) l=1:3
    pr.com=pc.svd.fdata(out) #
    out<-out[["data"]]
#   pr.com=prcomp(out)#,center=TRUE) #
    if (is.null(y)) stop("NO RESPONSE (y) INTRODUCED")
    if (nrow(out) != (length(y))) stop("ERROR IN THE DATA DIMENSIONS")
   }
 le=length(l)
 cor.y.pc=rep(NA,le)
 xxx=cbind(y,pr.com$x)
 cor.y.pc=round(cor(xxx[,c(1,l+1)]),3)[1,-1]
 names(cor.y.pc)=paste("cor(y,PC",l,")",sep="")
 if (draw){
  C<-match.call()
 lenC=length(C)
  j=1
  while (j<=lenC) {
     if (names(C)[j]=="ask") {ask=C[[j]];j=lenC +1}
          else {j=j+1;ask=FALSE}
  }
  dev.new()
  if (ask) {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          }
     else   par(mfrow=c(ceiling(le/2),2))
    for (i in 1:le)   plot(pr.com$x[,l[i]],y,main=paste("cor=",
        round(cor.y.pc[i],3)),xlab=colnames(pr.com$x)[l[i]],ylab="y",...)
  }
 return(cor.y.pc)
}
