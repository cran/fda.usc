fregre.pc=function (fdataobj, y, l = 1:3,norm=TRUE){
    if (!is.fdata(fdataobj))    fdataobj = fdata(fdataobj)
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x); np <- ncol(x);lenl = length(l)
    if (n != (length(y)))   stop("ERROR IN THE DATA DIMENSIONS")
    C <- match.call()
    if (is.null(rownames(x)))        rownames(x) <- 1:n
    ycen = y - mean(y)
    pc<-pc.svd.fdata(fdataobj,norm)
    xcen<-pc$fdataobj.cen
    if (length(l) == 1)  {
                  vs <- pc$rotation$data[l,]
                  Z<-pc$x[,l]
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
    beta.est$data <- as.numeric(beta.est$data)
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
    svd.fdata=pc,pf = pf,coefficients=object.lm$coefficients,residuals = object.lm$residuals,
    df = df,r2=r2, sr2 = sr2,H=H,fdataobj = fdataobj,y = y, l = l, pc=pc,lm=object.lm,Z=Z)
    class(out) = "fregre.fd"
    return(out)
}


