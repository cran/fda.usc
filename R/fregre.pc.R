####################################################################
####################################################################
fregre.pc=function (fdataobj, y, l =NULL,rn=0,weights=rep(1,len=n),...){
if (class(fdataobj)=="fdata.comp") {
    pc<-fdataobj
    fdataobj<-pc$fdataobj
    if (is.null(l))    {
       l<-pc$l
       }
    else if (length(l)>nrow(pc$rotation)) stop("Incorrect value for  argument l")
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]                                                    
   }
else {
 if (is.null(l)) l<- 1:3
 if (!is.fdata(fdataobj))    fdataobj = fdata(fdataobj)
#  omit<-omit.fdata(fdataobj,y)
#  fdataobj<-omit[[1]]
#  y<-omit[[2]]
  tt<-fdataobj[["argvals"]]
  x<-fdataobj[["data"]]
  pc<-fdata2pc(fdataobj,ncomp=max(l))
}

    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x); np <- ncol(x);lenl = length(l)
    if (is.null(rownames(x)))        rownames(x) <- 1:n
    X <-xcen<-pc$fdataobj.cen
    if (n != (length(y)))   stop("ERROR IN THE DATA DIMENSIONS")
    C <- match.call()
    ycen = y - mean(y)
    if (length(l) == 1)     vs<-pc$rotation$data
    else                    vs<-t(pc$rotation$data)
    scores<-Z<-(pc$x[,l])
    cnames<-colnames(pc$x)[l]
    df<-lenl+1
    J<-min(np,lenl)
    ymean<-mean(y)
    ycen<- y - ymean
    W<-diag(weights)    
 if (rn>0) {
    xmean<-pc$mean
    d<-pc$newd[l]
    D<-diag(d)
    diagJ<-diag(J)
    lenrn<-length(rn)
    scores<-cbind(rep(1,n),pc$x[,l])
    mat<-rn*diag(J+1)
    mat[1,1]<-0
#     xx2<-t(scores)%*%scores
     S<-solve(t(scores)%*%W%*%scores+mat)       
     Cinv<-S%*%t(scores)%*%W                    
     #comprobar y ver el siguiente
     coefs<-Cinv%*%y
     yp<-drop(scores%*%coefs)
     H<-scores%*%Cinv
     df<-traza(H)
     coefs<-drop(coefs)
    names(coefs)<-c("Intercept",cnames)
    beta.est<-coefs[-1]*pc$rotation[l]
#    beta.est$data<-apply(beta.est$data,2,sum)
     beta.est$data<-colSums(beta.est$data)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
    e<-y-yp
    rdf<-n-df
    sr2 <- sum(e^2)/ rdf
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
    GCV <- sum(e^2)/(n - df)^2
        object.lm = list()
        object.lm$coefficients <- coefs
        object.lm$residuals <- drop(e)
        object.lm$fitted.values <- yp
        object.lm$y <- y
        object.lm$rank <- df
        object.lm$df.residual <-  rdf
        Z=cbind(rep(1,len=n),Z)
        colnames(Z)[1] = "(Intercept)"
#        vcov2 = sr2 * Cinv
#        std.error = sqrt(diag(vcov2))
        std.error = sqrt(diag(S) *sr2)
        Vp<-sr2*S 
        t.value = coefs/std.error
        p.value = 2 * pt(abs(t.value), n - df, lower.tail = FALSE)
        coefficients <- cbind(coefs, std.error, t.value, p.value)
        colnames(coefficients) <- c("Estimate", "Std. Error",
            "t value", "Pr(>|t|)")
        class(object.lm) <- "lm"
#        b.est = b.est[-1]
#        names(b.est) <- rownames(coefficients)[-1]
    out <- list(call = C, beta.est = beta.est,coefficients=coefs,
    fitted.values =yp,residuals = e,H=H,df = df,r2=r2,GCV=GCV,
    sr2 = sr2,Vp=Vp,l = l,rn=rn,fdata.comp=pc,lm=object.lm,
    coefs=coefficients,fdataobj = fdataobj,y = y)
##################################
}
else {
#print("no rn")
    response = "y"
    dataf<-data.frame(y,Z,weights)
    colnames(dataf)<-c("y",cnames,"weights")
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf,"+",cnames[i],sep="")
    object.lm = lm(formula = pf,data=data.frame(dataf),weights=weights,x = TRUE, y = TRUE,...)
    beta.est<-object.lm$coefficients[2:(lenl+1)]*pc$rotation[l]
#    beta.est$data<-apply(beta.est$data,2,sum)
    beta.est$data<-colSums(beta.est$data)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
##    H<-diag(hat(Z, intercept = TRUE),ncol=n)
##    H<-lm.influence(object.lm, do.coef = FALSE)$hat o bien
##    I <- diag(1/(n*pc$lambdas[l]), ncol = lenl) #1/n

    Z=cbind(rep(1,len=n),Z)
    S=solve(t(Z)%*%W%*%Z)
    H<-Z%*%S%*%t(Z)
    e<-object.lm$residuals
#    H<-diag(hat(Z, intercept = TRUE),ncol=n)
#     H<-scores%*%solve(t(scores)%*%scores+mat)%*%t(scores)
#     df<-traza(H)
     df<-n- object.lm$df
     sr2 <- sum(e^2)/(n - df)
     Vp<-sr2*S 
     r2 <- 1 - sum(e^2)/sum(ycen^2)
     r2.adj<- 1 - (1 - r2) * ((n -    1)/(n-df))
     GCV <- sum(e^2)/(n - df)^2
  out <- list(call = C, beta.est = beta.est,coefficients=object.lm$coefficients,
  fitted.values =object.lm$fitted.values,residuals = e,H=H,df = df,r2=r2,GCV=GCV,
  sr2 = sr2,Vp=Vp, l = l,rn=rn,fdata.comp=pc,lm=object.lm,weights=weights,XX=Z,
  fdataobj = fdataobj,y = y)
  }
 class(out) = "fregre.fd"
 return(out)
}
####################################################################
####################################################################

