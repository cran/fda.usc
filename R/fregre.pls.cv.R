fregre.pls.cv=function (fdataobj, y, kmax=8, criteria = "SIC",...) {
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
    cv.opt1 = Inf;    pc.opt1 = NA
    ind =1:kmax
    num.pc = nrow(c)
    l = l2 = list()
    ck = 1
    max.c = length(c)
    tab = list("AIC", "AICc","SIC", "SICc","CV")
    type.i = pmatch(criteria, tab)
    if (is.na(type.i))     stop("Error: incorrect criteria")
    else {
    if (type.i < 5) {
        cv.AIC <- rep(NA, kmax)
        for (j in 1:kmax) {
            out = fregre.pls(fdataobj, y, l = 1:ind[j],...)
            s2 <- sum(out$residuals^2)/n  #(n-ck)
            ck=ind[j]+1
            if (criteria == "AIC") {
                cv.AIC[j] <- log(s2) + 2 * (ck+1)/n
            }
            else if (criteria == "AICc") {
                cv.AIC[j] <- log(s2) + 2 * (ck+1)/(n - ck - 2)
            }
            else if (criteria == "SIC") {
                cv.AIC[j] <- log(s2) + log(n) * ck/n
                }
            else if (criteria == "SICc") {
                cv.AIC[j] <- log(s2) + log(n) * ck/(n-ck-2)
            }
    }
    min.AIC = min(cv.AIC)
    pc.opt <- 1:ind[which.min(cv.AIC)]
    }
####
    else {
    pb=txtProgressBar(min=0,max=kmax,width=50,style=3)
    cv.AIC <- rep(NA, kmax)
    for (j in 1:kmax) {
        setTxtProgressBar(pb,j-0.5)
        residuals =residuals2<- rep(NA, n)
        ck=ind[j]
        for (i in 1:n){
            out = fregre.pls(fdataobj[-i,], y[-i],l = 1:ind[j],...)
            residuals[i] <- y[i] - predict(out,fdataobj[i,])
        }
        cv.AIC[j] <- sum(residuals^2)/(n-ck)
        min.AIC = min(cv.AIC)
        pc.opt <- 1:ind[which.min(cv.AIC)]
    }
    close(pb)
    }
###
    }
    names(cv.AIC) = paste("PLS",1:kmax , sep = "")
    fregre=fregre.pls(fdataobj,y,l=pc.opt,...)
    return(list("fregre.pls"=fregre,pls.opt = pc.opt, MSC.min = min.AIC,MSC = cv.AIC))
}

