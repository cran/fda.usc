#################################################################
#################################################################
fregre.pls.cv=function (fdataobj, y, kmax=8, criteria = "SICc",...) {
if (class(fdataobj)=="fdata.comp") {
    pc<-fdataobj
    fdataobj<-pc$fdataobj
    kmax<-nrow(pc$basis)
   }
else {
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
 omit<-omit.fdata(fdataobj,y)
 fdataobj<-omit[[1]]
 y<-omit[[2]]
 pc<-fdata2pls(fdataobj,y,kmax,...)
}
    x<-fdataobj[["data"]]
    tt<-fdataobj[["argvals"]]
    rtt<-fdataobj[["rangeval"]]
    n <- nrow(x);    nc <- ncol(x)
## if (is.logical(rn[1])) {
#   val<-log(.25*(pc$d[1]^2),base=2)
#   rn<-c(0,2^seq(0,val,len=10))
#   }
 # NO TENEMOS LA D PERO SI LOS SCORES
 # VAR(z)=d^2/n # calcular 1er autovalor lambda
     ind =1:kmax
    l = l2 = list()
    ck = 1
    tab = list("AIC", "AICc","SIC", "SICc","HQIC","rho","CV")
    type.i = pmatch(criteria, tab)
    pc2<-pc
    MSC.min<-Inf
    cv.AIC <- rep(NA,kmax)
    if (is.na(type.i))     stop("Error: incorrect criteria")
    else {
    if (type.i < 7) {
#        cv.AIC <- rep(NA, kmax)
        for (j in 1:kmax) {
            pc2$rotation<-pc$rotation[1:j]
            out = fregre.pls(pc2,y,...)
            ck<-out$df
            s2 <- sum(out$residuals^2)/n  #(n-ck)
            cv.AIC[j]<-switch(criteria,
              "AIC"=log(s2) + 2 * (ck)/n,
              "AICc"=log(s2) + 2 * (ck)/(n - ck - 2),
              "SIC"=log(s2) + log(n) * ck/n,
              "SICc"=log(s2) + log(n) * ck/(n-ck-2),
              "HQIC"=log(s2) + 2*log(log(n)) * ck/n,
              "rho"={A<-out$residuals;B<-1-diag(out$H)/n; D1<-(A/B)^2;sum(D1)})
     if ( MSC.min>cv.AIC[j]) {
       pc.opt<-j
       MSC.min= cv.AIC[j]
       }

#    min.AIC = min(cv.AIC)
#    pc.opt <- which.min(cv.AIC)
    }
    }
# CV criteria
    else {
         pc<-list()
         for (i in 1:n){
         xi<-fdataobj[-i]
         yi<-y[-i]
         pc[[i]] = fdata2pls(xi, yi,ncomp=kmax,...)
         }
        pc2<-pc
        for (j in 1:kmax) {
         residuals2<-rep(NA,n)
          for (i in 1:n){
            pc2[[i]]$rotation<-pc[[i]]$rotation[1:j]
            out = fregre.pls(pc2[[i]],y[-i],...)
            ck<-out$df
            a1<-out$coefficients[1]
            out$beta.est$data<-matrix(out$beta.est$data,nrow=1)
            b1<-inprod.fdata(fdata.cen(fdataobj[i],out$fdata.comp$mean)[[1]],out$beta.est)
            yp<- a1+b1
#            residuals[i] <- y[i] - yp
            residuals2[i] <- ((y[i] - yp)/(n-ck))^2
            }
#         cv.AIC[j] <- mean(residuals^2)/(n-j)^2###
         cv.AIC[j] <-sum(residuals2)/n

     if ( MSC.min>cv.AIC[j]) {
       pc.opt<-j
       MSC.min= cv.AIC[j]
       }
}     }   }
names(cv.AIC) = paste("PLS",1:kmax , sep = "")
#    rownames(cv.AIC) = paste("rn=",rn , sep = "")
    pc2$basis<-pc$rotation[1:pc.opt]
    fregre=fregre.pls(fdataobj,y,l=1:pc.opt,...)
    MSC.min = cv.AIC[pc.opt]
    return(list("fregre.pls"=fregre,pls.opt = 1:pc.opt,
    MSC.min = MSC.min,MSC = cv.AIC))
}
#################################################################
#################################################################

