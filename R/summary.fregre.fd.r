summary.fregre.fd=function(object,times.influ=3,times.sigma=3,draw=TRUE,...){
    x<-object$fdataobj$data
    t=object$fdataobj$argvals
    y<-object$y
    n=nrow(x)
 		up=mean(object$residuals)+times.sigma*sqrt(object$sr2)
    lo=mean(object$residuals)-times.sigma*sqrt(object$sr2)
    i.atypical=which(object$residuals>up|object$residuals<lo)
    lim.influ=traza(object$H)/n
    influence=diag(object$H)
    i.influence=which(influence>times.influ*lim.influ)
    if (length(i.influence) == 0) i.influence=NA
    if (length(i.atypical) == 0) i.atypical=NA
    if (object$call[[1]]=="fregre.pc") {
     cat(" *** Summary Functional Data Regression with Principal Components *** \n\n")
     object$lm$call<-object$call
     print(summary(object$lm))
     pr.x= apply(object$pc$x, 2, var)/sum(var(object$pc$x))
 cat("\n-With",length(object$l),"Principal Components is  explained ",round(sum(pr.x[object$l])*100
 ,2),"%\n of the variability of explicative variables. \n
-Variability for each  principal components -PC- (%):\n")
    print(round(pr.x[object$l] * 100, 2))
    }
    if (object$call[[1]]=="fregre.basis") {
     cat(" *** Summary Functional Data Regression with representation in Basis*** \n\n")
     if (object$lambda==0)     {object$lm$call<-object$call
                                print(summary(object$lm))}
      else  {
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefficients)
            cat("\n-R squared: ",object$r2)
            cat("\n-Residual variance: ",object$sr2,"\n")
            }
    }
  if (object$call[[1]]=="fregre.basis.cv") {
     cat(" *** Summary Functional Data Regression with representation in Basis*** \n\n")
      cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefficients)
            cat("\n-R squared: ",object$r2)
            cat("\n-Residual variance: ",object$sr2,"\n")
      cat("-Optimal Beta Basis: \n")
      print(object$basis.b.opt)
      cat("\n-Optimal lambda penalty=",object$lambda.opt,"\n")
#      print(object$Lfdobj)
    }
    if (object$call[[1]]=="fregre.np") {
     cat(" *** Summary  non-parametric Regression for Functional data*** \n\n")
     cat("-Call: ");    print(object$call)
     cat("\n-Bandwidth (h): ",object$h.opt)
    cat("\n-R squared: ",object$r2)
    cat("\n-Residual variance: ",object$sr2)

    }
  if (object$call[[1]]=="fregre.np.cv") {
     cat(" *** Summary  non-parametric Regression for Functional data*** \n\n")
     cat("-Call: ");    print(object$call)
     cat("\n-Bandwidth (h): ",object$h.opt)
    cat("\n-R squared: ",object$r2)
    cat("\n-Residual variance: ",object$sr2)

    }
     if (object$call[[1]]=="fregre.plm") {
     cat(" *** Summary semi-functional partial linear regression *** \n\n")
     cat("-Call: ");    print(object$call)
     cat("\n-Coefficients:  non functional covariates\n")
     print(object$coefficients)
      cat("\n-Bandwidth (h): ",object$h.opt)
    cat("\n-R squared: ",object$r2)
    cat("\n-Residual variance: ",object$sr2)

     }
    cat("\n-Names of possible atypical curves: ");
    if (is.na(i.atypical[1]))     cat("No atypical curves \n")
    else   if (length(i.atypical)<11)  cat(rownames(x)[i.atypical],"\n")
           else cat(rownames(x)[i.atypical[1:10]],
           "\n It prints only the 10 most atypical curves. \n")
    cat("-Names of possible influence curves : ");
    if (is.na(i.influence[1])) cat("No influence curves \n\n")
    else  if (length(i.influence)<11) cat(rownames(x)[i.influence],"\n\n")
     else cat(rownames(x)[i.influence[1:10]],
     "\n It prints only the 10 most influence curves \n\n")
    if (draw) {
      C<-match.call()
      lenC=length(C)
      j=1
      while (j<=lenC) {
        if (names(C)[j]=="ask") {
           ask=C[[j]]
           j=lenC +1             }
        else {      j=j+1
                    ask=FALSE             }
       }
       if (ask) {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          }
       else   par(mfrow=c(2,3))
 		plot(object$fitted.values,y,xlab="Fitted values",main=paste("R-squared=",
     round(object$r2,2)))
 		plot(object$fitted.values,object$residuals,ylab="Residuals",xlab="Fitted values",main="Residuals vs fitted.values")
    text(object$fitted.values[i.atypical],object$residuals[i.atypical],rownames(x)[i.atypical],cex=0.7)
    abline(h=mean(object$residuals),lwd=1,lty=2)
    abline(h=up,col=2,lwd=2,lty=2)
    abline(h=lo,col=2,lwd=2,lty=2)
#############
resid.sd=sqrt(abs(object$residuals/sd(object$residuals)))
main= "Scale-Location"
ylab23<-"Standardized residuals"
ylim <- c(0, max(resid.sd, na.rm = TRUE))
yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
plot(object$fitted.values,resid.sd, xlab = "Fitted values", ylab = yl, main = main,ylim = ylim)
 text(object$fitted.values[i.atypical],resid.sd[i.atypical],rownames(x)[i.atypical],cex=0.7)
#############
 		plot(diag(object$H),1:nrow(x),xlab="Leverage",ylab="Index.curves",
    main="Leverage")
text(diag(object$H)[i.influence],i.influence,rownames(x)[i.influence],cex=0.7)
 		abline(v=times.influ*lim.influ,col=2,lwd=2,lty=2)
#  	plot(density(object$residuals),main="Residuals")
    qqnorm(object$residuals,main="Residuals")
    boxplot(object$residuals,main="Residuals")
    }
    cat("\n")
return(invisible(list("Influence"=influence,"i.influence"=i.influence,"i.atypical"=i.atypical)))
}
