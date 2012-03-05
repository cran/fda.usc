fregre.pls.cv=function (fdataobj, y, kmax=8, criteria = "CV",...) {
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
    n <- nrow(x);    nc <- ncol(x)
#    cv.opt1 = Inf;    pc.opt1 = NA
    ind =1:kmax
#    num.pc = nrow(c)
    l = l2 = list()
    ck = 1
#    max.c = length(c)
    tab = list("AIC", "AICc","SIC", "SICc","HQIC","rho","CV")
    type.i = pmatch(criteria, tab)
    if (is.na(type.i))     stop("Error: incorrect criteria")
    else {
    if (type.i < 7) {
        cv.AIC <- rep(NA, kmax)
#        R<-correl(fdataobj)
        for (j in 1:kmax) {
            out = fregre.pls(fdataobj, y, l = 1:ind[j],...)
#             l2<-out$pls.fdata$res.pls$loadings
#             l3<-crossprod(l2,l2)
#           ck=traza((l3))
#print(ck)
#            ck=traza(t(l3)%*%l3)
             ck<-length(l)+1
            s2 <- sum(out$residuals^2)/n  #(n-ck)
            ck=ind[j]+1
            if (criteria == "AIC") {
                cv.AIC[j] <- log(s2) + 2 * (ck)/n
            }
            else if (criteria == "AICc") {
                cv.AIC[j] <- log(s2) + 2 * (ck)/(n - ck - 2)
            }
            else if (criteria == "SIC") {
                cv.AIC[j] <- log(s2) + log(n) * ck/n
                }
            else if (criteria == "SICc") {
                cv.AIC[j] <- log(s2) + log(n) * ck/(n-ck-2)
            }
            else if (criteria == "HQIC") {
                cv.AIC[j] <- log(s2) + 2*log(log(n)) * ck/n
            }
            else if (criteria == "rho") {
#              cv.AIC[j] <-   (traza(abs(out$residuals)*abs(R)*abs(out$residuals)))/(n-ck)
#              cv.AIC[j] <-   log(sum(out$residuals*abs(R)*out$residuals))/n +2*log(log(n)) * ck/n        }
#cv.AIC[j] <-   log(matrix(out$residuals,nrow=1)%*%abs(R)%*%out$residuals)/n +2*log(log(n)) * ck/n
#va bien si usamos las 3 siguientes lineas
#A<-out$residuals*abs(R)*out$residuals
#B<-1-diag(out$H);D1<-A/B
#cv.AIC[j] <-   log(sum(D1)/n) + log(n) * ck/(n-ck-2)


#    ss<-drop(y)*(1-out$H)*drop(y)
#       print(dim(ss))
#    s1<-solve(ss)
#   print(dim(s1))
A<-out$residuals
B<-1-diag(out$H)/n
D1<-(A/B)^2
cv.AIC[j] <-   sum(D1)#+ 2*traza(out$H)*out$sr2/n
                }
    }
    min.AIC = min(cv.AIC)
    pc.opt <- 1:ind[which.min(cv.AIC)]
    }
####
    else {
#    pb=txtProgressBar(min=0,max=kmax,width=50,style=3)
    cv.AIC <- rep(NA, kmax)
#    for (j in 1:kmax) {
#        setTxtProgressBar(pb,j-0.5)
#        residuals =residuals2<- rep(NA, n)
#       ck<-ind[j]
#        for (i in 1:n){
#            out = fregre.pls(fdataobj[-i,], y[-i],l = 1:ind[j],...)
#           residuals[i] <- y[i] - predict(out,fdataobj[i,])
#        }
#        cv.AIC[j] <- sum(residuals^2)/(n-ck)
        out = fregre.pls(fdataobj,y,l=1:kmax,...)
        cv.AIC <- out$fdata.comp$res.pls$validation$PRESS
        min.AIC = min(cv.AIC)
        pc.opt <- 1:which.min(cv.AIC)
#    }
#    close(pb)
#    }
    }   }
    names(cv.AIC) = paste("PLS",1:kmax , sep = "")
    fregre=fregre.pls(fdataobj,y,l=pc.opt,...)
    return(list("fregre.pls"=fregre,pls.opt = pc.opt, MSC.min = min.AIC,MSC = cv.AIC))
}

