fdata.cen=function(fdataobj){
if (!is.fdata(fdataobj))  fdataobj=fdata(fdataobj)
data<-fdataobj[["data"]]
 meanX <-func.mean(fdataobj)
 Xcen <- sweep(data,2,meanX[["data"]],FUN="-")
 Xcen<-fdata(Xcen,fdataobj[["argvals"]])
 return(list("Xcen"=Xcen,"meanX"=meanX))
}

#fdata.cen22 <- function(fdata){
#n <- dim(fdata)[1]
#unos <- rep(1,n)
#meanX <- (1/n) * t(fdata) %*% unos
#Xcen <- fdata - unos %*% t(meanX)
#return(list("Xcen"=Xcen,"meanX"=meanX))
#}




