semimetric.basis=function(DATA1, DATA2 = DATA1,nderiv=0,type.basis1=NULL,
nbasis1=NULL,type.basis2=NULL,nbasis2=NULL,...) {

 if (any(class(DATA1)=="fd")) {
   print("fd class")
   r=DATA1$basis[[3]]
   tt=seq(r[1],r[2],len=length(DATA1$fdnames$time))
   df1=deriv.fd(DATA1,nderiv)
   df2=deriv.fd(DATA2,nderiv)
   fd1=t(eval.fd(tt,df1))
   fd2=t(eval.fd(tt,df2))
   mdist=metric.lp(fd1,fd2,...)
  }
 else {
 if (!is.fdata(DATA1)) DATA1=fdata(DATA1)
 if (!is.fdata(DATA2)) DATA2=fdata(DATA2)
  tt=DATA1[["argvals"]]
  rtt<-DATA1[["rangeval"]]
  DATA1<-DATA1[["data"]]
  DATA2<-DATA2[["data"]]
#   print("Raw class")
   np=ncol(DATA1)
   if (is.null(nbasis1)) {
#       nbasis1=floor(np/6)
       nbasis1=ifelse(floor(np/3) > floor((np - nderiv - 4)/2),
       floor((np - nderiv - 4)/2), floor(np/3))
       }
   if (is.null(nbasis2)) nbasis2=nbasis1
   as <- list()
   bs <- list()
   as[[1]] <- rtt
   bs[[1]] <- rtt
   names(as)[[1]]<-"rangeval"
   names(bs)[[1]]<-"rangeval"
   as[[2]] <- nbasis1
   names(as)[[2]]<-"nbasis"
   C <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m<-match(c("DATA1", "DATA2","nderiv","type.basis1","nbasis1","type.basis2","nbasis2"),names(mf),0L)
   imetric <- m[4]
   imetric2 <- m[6]
   if (imetric == 0) {
        type.basis1="bspline"
        a1 <- create.bspline.basis
        len.metric <- length(formals(a1))
        vv <- array(0, dim = c(len.metric))    }
   else {  a1 <- paste('create.',type.basis1,'.basis',sep="")
        len.metric <- length(formals(a1))
        vv <- array(0, dim = c(len.metric)) }
  ii <- imetric + 1
  if (C[ii] != "NULL()") {
        ind.m <- 3
        while (C[ii] != "NULL()" && ind.m <= len.metric) {
            aa <- any(names(C) == names(formals(a1))[ind.m])
            if (aa) {
                vv[ind.m] <- which(names(C) == names(formals(a1)[ind.m]))
                ii <- ii + 1
                as[[ind.m]] <- C[[vv[ind.m]]]
                names(as)[[ind.m]]<-names(formals(a1)[ind.m])            }
#            else {                 as[[ind.m]] <- formals(a1)[[ind.m]]   }
            ind.m <- ind.m + 1            }
  }
 b1.1<- do.call(a1, as)
 if (imetric2 == 0) {
         b1 <-  paste('create.',type.basis1,'.basis',sep="")
        len.metric <- length(formals(b1))
        vv <- array(0, dim = c(len.metric))
}
else {  b1 <- paste('create.',type.basis2,'.basis',sep="")
        len.metric <- length(formals(b1))
        vv <- array(0, dim = c(len.metric)) }
 ii <- imetric2 + 1
 if (C[ii] != "NULL()") {
        ind.m <- 3
        while (C[ii] != "NULL()" && ind.m <= len.metric) {
            aa <- any(names(C) == names(formals(b1))[ind.m])
            if (aa) {
                vv[ind.m] <- which(names(C) == names(formals(b1)[ind.m]))
                ii <- ii + 1
                bs[[ind.m]] <- C[[vv[ind.m]]]
                names(bs)[[ind.m]]<-names(formals(b1)[ind.m])            }
            else { as[[ind.m]] <- formals(b1)[[ind.m]]}
            ind.m <- ind.m + 1}
  }

   bs[[2]] <- nbasis2
   names(bs)[[2]]<-'nbasis'
   b1.2<- do.call(b1, bs)
   class(DATA1)="matrix"
   class(DATA2)="matrix"
   fd1.1 <- Data2fd(argvals=tt,y=t(DATA1),basisobj=b1.1)
   fd1.2 <- Data2fd(argvals=tt,y=t(DATA2),basisobj=b1.2)
   r=range(tt)
   df1=deriv.fd(fd1.1,nderiv)
   df2=deriv.fd(fd1.2,nderiv)
   fd1=t(eval.fd(tt,df1))
   fd2=t(eval.fd(tt,df2))
   mdist=metric.lp(fd1,fd2,...)
}
mdist
}




