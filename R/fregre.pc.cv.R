fregre.pc.cv=function (fdataobj, y, kmax=8, criteria = "SIC",...)
{
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
    n <- nrow(x)
    nc <- ncol(x)
    cv.opt1 = Inf
    pc.opt1 = NA
    c = matrix(1:kmax, nrow = 1)
    num.pc = nrow(c)
    l = l2 = list()
    ck = 1
    max.c = length(c)
    c0 = 1:kmax
    use = rep(FALSE, kmax)
    tab = list("AIC", "AICc","SIC", "SICc","CV")
    type.i = pmatch(criteria, tab)
    if (is.na(type.i))     stop("Error: incorrect criteria")
    else {
    if (type.i < 5) {
    for (k in 1:kmax) {
        cv.AIC <- rep(NA, max.c)
        for (j in 1:max.c) {
            out = fregre.pc(fdataobj, y, l = c[, j],...)
            s2 <- sum(out$residuals^2)/n  #(n-ck)
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
        pc.opt1 <- c[, which.min(cv.AIC)]
        l[[k]] = pc.opt1[k]
        l2[[k]] = min.AIC
        use[pc.opt1[k]] = TRUE
        l[[k + 1]] = c0[use == FALSE]
        c = t(expand.grid(l))
        ck=nrow(c)+1
        max.c = ncol(c)
    }
    }
    else {
    pb=txtProgressBar(min=0,max=kmax,width=50,style=3)
    for (k in 1:kmax) {
      setTxtProgressBar(pb,k-0.5)
      cv.AIC <- rep(NA, max.c)
      for (j in 1:max.c) {
       residuals <- rep(NA, n)
        for (i in 1:n){
            print("no llega")
           out = fregre.pc(fdataobj[-i,], y[-i], l = c[, j])
           y.est=out$a.est*rep(1,n)+x%*%out$beta.est[["data"]][1,]/(nc-1)      ####
           residuals[i] <- y[i] - y.est[i]
          }
         cv.AIC[j] <- sum(residuals^2)/(n-ck)
        }
        min.AIC = min(cv.AIC)
        pc.opt1 <- c[, which.min(cv.AIC)]
        l[[k]] = pc.opt1[k]
        l2[[k]] = min.AIC
        use[pc.opt1[k]] = TRUE
        l[[k + 1]] = c0[use == FALSE]
        c = t(expand.grid(l))
        ck=nrow(c)+1
        max.c = ncol(c)
         setTxtProgressBar(pb,k)
    }
    close(pb)
    }
    mn = which.min(l2)
    pc.opt = pc.opt1[1:mn]
    MSC = t(l2)
    }
    for (i in 1:kmax) names(MSC)[i] = paste("PC", t(l[i]), sep = "")
    fregre=fregre.pc(fdataobj,y,l=pc.opt,...)
    return(list("fregre.pc"=fregre,pc.opt = pc.opt, MSC.min = l2[[mn]], pc.order =
t(l[1:kmax]),MSC = MSC))
}
