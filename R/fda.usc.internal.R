"[.fdata"<-function(fdataobj, i = TRUE, j = TRUE,drop=FALSE) {
if (is.numeric(j) && j==1 && length(j)==1)
fdataobj[["data"]]<-matrix(fdataobj[["data"]][i,j],nrow=1)
else  fdataobj[["data"]]<-fdataobj[["data"]][i,j,drop=drop]
fdataobj[["argvals"]]<-fdataobj[["argvals"]][j]
fdataobj
}
################################################################################

"==.fdata"<-function(fdata1,fdata2){
fdataequal<-TRUE
 if (!(all(fdata1[["data"]] == fdata2[["data"]]))) {
        fdataequal <- FALSE
        print("No equal data matrix")
    }
 if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]]))) {
        fdataequal <- FALSE
        print("No equal argvals vector")
    }
 return(fdataequal)
}
################################################################################

"+.fdata"<-function(fdata1,fdata2){
inhe.fdata1<- inherits(fdata1, "fdata")
inhe.fdata2<- inherits(fdata2, "fdata")
if (!inhe.fdata1 && !inhe.fdata2)
  stop("Neither argument for * is a functional data object fdata.")
if (inhe.fdata1 && inhe.fdata2) {
 if (!(all(dim(fdata1[["data"]])==dim(fdata1[["data"]])))) stop("Error in data dimenstion")
 if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]])))  stop("Error in argvals")
 fdataobj<-fdata1
 fdataobj[["data"]]<-fdata1[["data"]]+fdata2[["data"]]
}
if (!inhe.fdata1 && inhe.fdata2) {
 fdataobj<-fdata2
 fdataobj[["data"]]<-fdata1+fdata2[["data"]]
}
if (inhe.fdata1 && !inhe.fdata2) {
 fdataobj<-fdata1
 fdataobj[["data"]]<-fdata1[["data"]]+fdata2
}
fdataobj
}
################################################################################
"-.fdata"<-function(fdata1,fdata2){
inhe.fdata1<- inherits(fdata1, "fdata")
inhe.fdata2<- inherits(fdata2, "fdata")
if (!inhe.fdata1 && !inhe.fdata2)
  stop("Neither argument for * is a functional data object fdata.")
if (inhe.fdata1 && inhe.fdata2) {
 if (!(all(dim(fdata1[["data"]])==dim(fdata1[["data"]])))) stop("Error in data dimenstion")
 if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]])))  stop("Error in argvals")
 fdataobj<-fdata1
 fdataobj[["data"]]<-fdata1[["data"]]-fdata2[["data"]]
}
if (!inhe.fdata1 && inhe.fdata2) {
 fdataobj<-fdata2
 fdataobj[["data"]]<-fdata1-fdata2[["data"]]
}
if (inhe.fdata1 && !inhe.fdata2) {
 fdataobj<-fdata1
 fdataobj[["data"]]<-fdata1[["data"]]-fdata2
}
fdataobj
}

################################################################################

"*.fdata"<-function(fdata1,fdata2){
inhe.fdata1<- inherits(fdata1, "fdata")
inhe.fdata2<- inherits(fdata2, "fdata")
if (!inhe.fdata1 && !inhe.fdata2)
  stop("Neither argument for * is a functional data object fdata.")
if (inhe.fdata1 && !inhe.fdata2) {
 fdataobj<-fdata1
 fdataobj[["data"]]<-fdata1[["data"]]*fdata2
}
if (!inhe.fdata1 && inhe.fdata2) {
 fdataobj<-fdata2
 fdataobj[["data"]]<-fdata1*fdata2[["data"]]
}
if (inhe.fdata1 && inhe.fdata2) {
 if (!(all(dim(fdata1[["data"]])==dim(fdata1[["data"]])))) stop("Error in data dimenstion")
 if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]])))  stop("Error in argvals")
 fdataobj<-fdata1
 fdataobj[["data"]]<-fdata1[["data"]]*fdata2[["data"]]
}
fdataobj
}
################################################################################
"/.fdata"<-function(fdata1,fdata2){
inhe.fdata1<- inherits(fdata1, "fdata")
inhe.fdata2<- inherits(fdata2, "fdata")
if (!inhe.fdata1 && !inhe.fdata2)
  stop("Neither argument for * is a functional data object fdata.")
if (inhe.fdata1 && !inhe.fdata2) {
 fdataobj<-fdata1
 fdataobj[["data"]]<-fdata1[["data"]]/fdata2
}
if (!inhe.fdata1 && inhe.fdata2) {
 fdataobj<-fdata2
 fdataobj[["data"]]<-fdata1/fdata2[["data"]]
}
if (inhe.fdata1 && inhe.fdata2) {
 if (!(all(dim(fdata1[["data"]])==dim(fdata1[["data"]])))) stop("Error in data dimenstion")
 if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]])))  stop("Error in argvals")
 fdataobj<-fdata1
 fdataobj[["data"]]<-fdata1[["data"]]/fdata2[["data"]]
}
fdataobj
}


################################################################################
dim.fdata<-function(x) {dim(x[["data"]])}
ncol.fdata<-function(fdataobj){ncol(fdataobj[["data"]])}
nrow.fdata<-function(fdataobj){nrow(fdataobj[["data"]])}
#"length.fdata"<-function(fdataobj){
#length(fdataobj[["argvals"]])}

#getS3method("[","fd")
#getS3method("==","basisfd")

#getS3method("[","fdata")
#a<-tecator[[1]][1,]

c.fdata<-function(...) {
    C=match.call()
    fdatalist <- list(...)
    n <- length(fdatalist)
    v<-rep(FALSE,len=n)
    fdata1 <- fdatalist[[1]]
    if (n == 1)  return(fdata1)
    if (is.vector(fdata1$data))  {
         fdata1$data=matrix(fdata1$data,nrow=1)
         v[1]<-TRUE
         }
    data<- fdata1$data
    dimdata <- dim(data)
    ndim <- length(dimdata)
    argvals <- fdata1$argvals
    rangeval <- fdata1$rangeval
    names <- fdata1$names
    if (!inherits(fdata1, "fdata")) stop("Objects must be of class fdata")
    for (j in (2:n)) {
        fdataj <- fdatalist[[j]]
        if (is.vector(fdataj$data))  {
           fdataj$data=matrix(fdataj$data,nrow=1)
           v[j]<-TRUE
           fdatalist[[j]]<-fdataj
           }
        if (!inherits(fdataj, "fdata"))
            stop("Objects must be of class fdata")
        if (any(unlist(fdataj$argvals) != unlist(argvals)))
            stop("Objects must all have the same argvals")
        if (any(unlist(fdataj$rangeval) != unlist(rangeval)))
            stop("Objects must all have the same rangeval")
        if (any(unlist(fdataj$names) != unlist(names)))        {
            print("Concatenate main names")
            names$main<-paste(names$main,"_",fdataj$names$main,sep="")
            }
        if (length(dim(fdataj$data)) != ndim)
            stop("Objects must all have the same number of multiple functions")
    }
    if (ndim == 2) {
        for (j in 2:n) {
           fdataj <- fdatalist[[j]]
           dataj <- fdataj$data
           dd<-C[[j+1]]
           if (v[j])    rownames(dataj)<-deparse(substitute(dd))
           data <- rbind(data, dataj)
        }
        dd1<-C[[2]]
       if (v[1])         rownames(data)[1]<-deparse(substitute(dd1))
    }
    concatfdata <- fdata(data, argvals,rangeval, names)
    return(concatfdata)
}

