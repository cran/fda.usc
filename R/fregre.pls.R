#################################################################
#################################################################
fregre.pls=function(fdataobj, y=NULL, l = NULL,...){
if (class(fdataobj)=="fdata.comp") {
    pc<-fdataobj
    fdataobj<-pc$fdataobj
    if  (is.null(l)) l<-1:nrow(pc$rotation)
    if  (is.null(y)) y<-pc$y
    else if (all(y!=pc$y)) warning("y is different from that calculated on the pls basis")
     lenl<-max(l)
   }
else {
 if (is.null(l)) l<- 1:3
# omit<-omit.fdata(fdataobj,y)
# fdataobj<-omit[[1]]
# y<-omit[[2]]
 lenl<-max(l)
 pc<-fdata2pls(fdataobj,y,lenl,...)
}
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x); J<-np <- ncol(x);#lenl = length(l)
    if (n != (length(y)))   stop("ERROR IN THE DATA DIMENSIONS")
    C <- match.call()
    if (is.null(rownames(x)))        rownames(x) <- 1:n
    ycen = y - mean(y)
    if (length(l) == 1)      vs <- pc$rotation$data
    else                     vs <- t(pc$rotation$data)
    Z<-(pc$x[,l])
    xcen<-pc$fdataobj.cen
    cnames<-colnames(pc$x)[l]
    response = "y"
    XX<-data.frame(y,Z)
    colnames(XX)<-c("y",cnames)
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf,"+",cnames[i],sep="")
#   pf=paste(pf,"-1")
    object.lm = lm(formula = pf, data =XX , ...)
    beta.est<-object.lm$coefficients[2:(lenl+1)]*pc$rotation[l]
#    beta.est$data<-apply(beta.est$data,2,sum)
     beta.est$data<-colSums(beta.est$data)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
#    pc$df
#    H<-diag(hat(Z, intercept = TRUE),ncol=n)
 # H2<-lm.influence(object.lm, do.coef = T)$hat# o bien
 #    I <- diag(1/(n*pc$lambdas[l]), ncol = lenl) #1/n
    Z=cbind(rep(1,len=n),Z)
    S=solve(t(Z)%*%Z)
    H<-Z%*%S%*%t(Z)
    e<-object.lm$residuals
#    df = traza(H)
    df<-pc$df[lenl]+1
    rdf<-n-df
    object.lm$df.residual<-rdf
    sr2 <- sum(e^2)/rdf
    Vp<-sr2*S
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
    GCV <- sum(e^2)/rdf^2
    GCV <- sum(e^2)/(n - df)^2
#        object.lm$coefficients <- coefs
#        object.lm$rank <- df
        Z=cbind(rep(1,len=n),Z)
        colnames(Z)[1] = "(Intercept)"
        std.error = sqrt(diag(S) *sr2)
        t.value =object.lm$coefficients/std.error
        p.value = 2 * pt(abs(t.value), n - df, lower.tail = FALSE)
        coefficients <- cbind(object.lm$coefficients, std.error, t.value, p.value)
        colnames(coefficients) <- c("Estimate", "Std. Error",
            "t value", "Pr(>|t|)")
        class(object.lm) <- "lm"
 out <- list(call = C, beta.est = beta.est,coefficients=object.lm$coefficients,
 fitted.values =object.lm$fitted.values, residuals = object.lm$residuals,
 H=H,df = df,r2=r2, GCV=GCV,sr2 = sr2,Vp=Vp, l = l, fdata.comp=pc, coefs=coefficients,
 lm=object.lm,fdataobj = fdataobj,y = y,XX=XX)
    class(out) = "fregre.fd"
    return(out)
}
#################################################################
#################################################################
