#################################################################
#################################################################
fregre.pc.cv=function (fdataobj, y, kmax=8,rn=0,criteria = "SIC",...) {
if (class(fdataobj)=="fdata.comp") {
    pc<-fdataobj
    fdataobj<-pc$fdataobj
    if (is.null(kmax))    {
       kmax<-max(pc$l)
       }
    else if (kmax>nrow(pc$rotation)) stop("Incorrect value for  argument l")
#    kmax<-length(pc$l)
    tt<-fdataobj[["argvals"]]
    x<-fdataobj[["data"]]
    np<-ncol(x)
   }
else {
  if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  tt<-fdataobj[["argvals"]]
  x<-fdataobj[["data"]]
  np<-ncol(x)
  pc<-fdata2pc(fdataobj,ncomp=kmax,...)
}
if (is.logical(rn[1])) {
   val<-log(.25*(pc$d[1]^2),base=2)
   rn<-c(0,2^seq(0,val,len=10))
   }
if (is.null(names(y))) names(y)<-1:length(y)
#x<-fdataobj[["data"]]
#tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
n <- nrow(x);    #nc <- ncol(x)
cv.opt1 = Inf;    pc.opt1 = NA
c1 = matrix(1:kmax, nrow = 1)
num.pc = nrow(c1)
l = l2 = list()
#ck = 1
max.c = length(c1)
c0 = 1:kmax
use = rep(FALSE, kmax)
tab = list("AIC", "AICc","SIC", "SICc","rho","CV")
type.i = pmatch(criteria, tab)
#pc2<-pc
lenrn<-length(rn)
MSC3<-list()
pc.opt2 <- matrix(NA,nrow=lenrn,ncol=kmax)
rownames(pc.opt2)<-paste("rn=",zapsmall(signif(rn)),sep="")
colnames(pc.opt2)<-paste("PC(",1:kmax,")",sep="")
MSC2<-pc.opt2
MSC.min<-Inf
min.rn<-rn[1]
if (is.na(type.i))     stop("Error: incorrect criteria")
else {
    if (type.i < 6) {
    for (r in 1:lenrn) {
    cv.opt1 = Inf
    pc.opt1 = NA
    l = l2 = list()
     c1 = matrix(1:kmax, nrow = 1)
         num.pc = nrow(c1)
     max.c = length(c1)
     c0 = 1:kmax
     use = rep(FALSE, kmax)
     pc2<-pc
     for (k in 1:kmax) {
                cv.AIC <- rep(NA, max.c)
                for (j in 1:max.c) {
                  pc2$rotation <- pc$rotation#[c1[, j]]
                  pc2$l <- pc$l[c1[, j]]
                  out = fregre.pc(pc2, y,l=c1[, j],rn=rn[r],...)
                  ck<-out$df
                  s2 <- sum(out$residuals^2)/n
                  cv.AIC[j]<-switch(criteria,
              "AIC"=log(s2) + 2 * (ck)/n,
              "AICc"=log(s2) + 2 * (ck)/(n - ck - 2),
              "SIC"=log(s2) + log(n) * ck/n,
              "SICc"=log(s2) + log(n) * ck/(n-ck-2),
              "HQIC"=log(s2) + 2*log(log(n)) * ck/n,
              "rho"={A<-out$residuals;B<-1-diag(out$H)/n; D1<-(A/B)^2;sum(D1)})
                }
                #peta en 2!!!
                min.AIC = min(cv.AIC)
                pc.opt1 <- c1[, which.min(cv.AIC)]
                l[[k]] = pc.opt1[k]
                l2[[k]] = min.AIC
                use[pc.opt1[k]] = TRUE
                l[[k + 1]] = c0[use == FALSE]
                c1 = t(expand.grid(l))
                ck = nrow(c1) + 1
                max.c = ncol(c1)
            }
     mn = which.min(l2)
     MSC = as.numeric(l2)
     if ( MSC.min>MSC[mn]) {
       min.rn<-r
       MSC.min = MSC[mn]
       pc.opt3<-pc.opt1[1:mn]
       }
    pc.opt = pc.opt1[1:mn]
    MSC2[r,]<-MSC
    pc.opt2[r,]<-pc.opt1
    }
    }
#### CV criteria
    else {
    pb=txtProgressBar(min=0,max=lenrn,width=50,style=3)
    pcl<-list()
    for (i in 1:n) {pcl[[i]]<-fdata2pc(fdataobj[-i,],ncomp=kmax,...)}
    for (r in 1:lenrn) {
     setTxtProgressBar(pb,r-0.5)
     cv.opt1 = Inf
     pc.opt1 = NA
     l = l2 = list()
     c1 = matrix(1:kmax, nrow = 1)
     num.pc = nrow(c1)
     max.c = length(c1)
     c0 = 1:kmax
     use = rep(FALSE, kmax)
       for (k in 1:kmax) {
      cv.AIC <- rep(NA, max.c)
      cv.AIC2 <- matrix(NA,nrow=max.c,ncol=lenrn)
      rownames(cv.AIC2)<-1:max.c
      for (j in 1:max.c) {
        residuals2<- rep(NA, n)
        maxk<-max(c1[, j])
          for (i in 1:n){
            pc2<-pcl[[i]]
            pc2$rotation<-pcl[[i]]$rotation#[c1[,j]]
            pc2$l<-pcl[[i]]$l[c1[,j]]
            out = fregre.pc(pc2,y[-i],l=c1[, j],rn=rn[r],...) #####
#            out = fregre.pc(fdataobj[-i,],y[-i],l=c1[,j],rn=rn[r],...)
            ck<-out$df
            residuals2[i] <- ((y[i] - predict(out,fdataobj[i,]))/(n-ck))^2
#            cat(ck)
            }
          cv.AIC[j] <-sum(residuals2)/n
          }
                min.AIC = min(cv.AIC)
                pc.opt1 <- c1[, which.min(cv.AIC)]
                l[[k]] = pc.opt1[k]
                l2[[k]] = min.AIC
                use[pc.opt1[k]] = TRUE
                l[[k + 1]] = c0[use == FALSE]
                c1 = t(expand.grid(l))
                ck = nrow(c1) + 1
                max.c = ncol(c1)
            }
     mn = which.min(l2)
     MSC = as.numeric(l2)
     if ( MSC.min>MSC[mn]) {
       min.rn<-r
       MSC.min = MSC[mn]
       pc.opt3<-pc.opt1[1:mn]
       }
    pc.order<-names(MSC)
    pc.opt = pc.opt1[1:mn]
    MSC2[r,]<-MSC
    pc.opt2[r,]<-pc.opt1
    setTxtProgressBar(pb,r)
    }
    close(pb)
    }
 mn = which.min(l2)
 MSC = as.numeric(l2)
 names(pc.opt3)<-paste("PC", pc.opt3, sep = "")
 rn.min<-rn[min.rn]
 fregre=fregre.pc(fdataobj,y,l=pc.opt,rn=rn.min,...)
 return(list("fregre.pc"=fregre,pc.opt = pc.opt3,rn.opt=rn.min,PC.order=pc.opt2,
 MSC.order=MSC2))
}
}
#################################################################
#################################################################

