#1. preg a manolo si max(l) o sequencial!!!
#2. no hace falta hacer object.lm<-lm(df....etc etc ) o asi solucionamos 1
fregre.pls=function (fdataobj, y, l = 1:3,...){
if (!is.fdata(fdataobj))    fdataobj = fdata(fdataobj)
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
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x); np <- ncol(x);lenl = length(l)
    if (n != (length(y)))   stop("ERROR IN THE DATA DIMENSIONS")
    C <- match.call()
    if (is.null(rownames(x)))        rownames(x) <- 1:n
    ycen = y - mean(y)
    pc<-pls.fdata(fdataobj,ycen,max(l),...)
    xcen<-pc$fdataobj.cen
     if (length(l) == 1)  {
                  vs <- pc$rotation$data[l,]
                  Z<-pc$x[,1:l]
                  }
    else {
                  vs <- t(pc$rotation$data[l,])
                  Z<-(pc$x[,l])
                  }
    cnames<-colnames(pc$x)[l]
    response = "y"
    df<-data.frame(y,Z)
    colnames(df)<-c("y",cnames)
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf,"+",cnames[i],sep="")
#   pf=paste(pf,"-1")
    object.lm = lm(formula = pf, data =df , x = TRUE,y = TRUE)
    beta.est<-object.lm$coefficients[2:(lenl+1)]*pc$rotation[l,]
    beta.est$data<-apply(beta.est$data,2,sum)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
#    H<-diag(hat(Z, intercept = TRUE),ncol=n)
# H<-lm.influence(object.lm, do.coef = FALSE)$hat o bien
#    I <- diag(1/(n*pc$lambdas[l]), ncol = lenl) #1/n
    Z=cbind(rep(1,len=n),Z)
    S=solve(t(Z)%*%Z)
    H<-Z%*%S%*%t(Z)
#    xmean=apply(x,2,mean)
    e<-object.lm$residuals
    df = lenl + 1
    sr2 <- sum(e^2)/(n - df)
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    out <- list(call = C, beta.est = beta.est, fitted.values =object.lm$fitted.values,
    pls.fdata=pc,coefficients=object.lm$coefficients,residuals = object.lm$residuals,
    df = df,r2=r2, sr2 = sr2,H=H,fdataobj = fdataobj,y = y, l = l, lm=object.lm)
    #,Z=Z,pf = pf)
    class(out) = "fregre.fd"
    return(out)
}

