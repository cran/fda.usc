################################################################################

"+.fdata" <- function(fdata1,fdata2){
  inhe.fdata1 <-  inherits(fdata1, "fdata")
  inhe.fdata2 <-  inherits(fdata2, "fdata")
  if (!inhe.fdata1 && !inhe.fdata2)
    stop("Neither argument for + is a functional data object fdata.")
  if (inhe.fdata1 && inhe.fdata2) {
    if (ncol(fdata1)!=ncol(fdata2)) stop("Error in columns dimensions")
    if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]])))  stop("Error in argvals")
    n1=nrow(fdata1)
    n2=nrow(fdata2)
    if (n1==n2){
      fdataobj <- fdata1
      fdataobj[["data"]] <- fdata1[["data"]]+fdata2[["data"]] } 
    else {
      if (n1<n2) {
        if (n2%%n1!=0) stop("Number of rows not compatible")
        fdataobj=fdata2
        fdataobj[["data"]]=fdata2[["data"]]+matrix(rep(t(fdata1$data),n2%/%n1),ncol=ncol(fdata1),byrow=TRUE)
      } else {
        if (n1%%n2!=0) stop("Number of rows not compatible")
        fdataobj=fdata1
        fdataobj[["data"]]=fdata1[["data"]]+matrix(rep(t(fdata2$data),n1%/%n2),ncol=ncol(fdata2),byrow=TRUE)
      }
    }
  }
  if (!inhe.fdata1 && inhe.fdata2) {
    fdataobj <- fdata2
    fdataobj[["data"]] <- fdata1+fdata2[["data"]]
  }
  if (inhe.fdata1 && !inhe.fdata2) {
    fdataobj <- fdata1
    fdataobj[["data"]] <- fdata1[["data"]]+fdata2
  }
  fdataobj
}
################################################################################

"-.fdata" <- function(fdata1,fdata2){
  inhe.fdata1 <-  inherits(fdata1, "fdata")
  inhe.fdata2 <-  inherits(fdata2, "fdata")
  if (!inhe.fdata1 && !inhe.fdata2)
    stop("Neither argument for - is a functional data object fdata.")
  if (inhe.fdata1 && inhe.fdata2) {
    if (ncol(fdata1)!=ncol(fdata2)) stop("Error in columns dimensions")
    if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]])))  stop("Error in argvals")
    n1=nrow(fdata1)
    n2=nrow(fdata2)
    if (n1==n2){
      fdataobj <- fdata1
      fdataobj[["data"]] <- fdata1[["data"]]-fdata2[["data"]] } 
    else {
      if (n1<n2) {
        if (n2%%n1!=0) stop("Number of rows not compatible")
        fdataobj=fdata2
        fdataobj[["data"]]=matrix(rep(t(fdata1$data),n2%/%n1),ncol=ncol(fdata1),byrow=TRUE)-fdata2[["data"]]
      } else {
        if (n1%%n2!=0) stop("Number of rows not compatible")
        fdataobj=fdata1
        fdataobj[["data"]]=fdata1[["data"]]-matrix(rep(t(fdata2$data),n1%/%n2),ncol=ncol(fdata2),byrow=TRUE)
      }
    }
  }
  if (!inhe.fdata1 && inhe.fdata2) {
    fdataobj <- fdata2
    fdataobj[["data"]] <- fdata1-fdata2[["data"]]
  }
  if (inhe.fdata1 && !inhe.fdata2) {
    fdataobj <- fdata1
    fdataobj[["data"]] <- fdata1[["data"]]-fdata2
  }
  fdataobj
}
################################################################################

"*.fdata" <- function(fdata1,fdata2){
  inhe.fdata1 <-  inherits(fdata1, "fdata")
  inhe.fdata2 <-  inherits(fdata2, "fdata")
  if (!inhe.fdata1 && !inhe.fdata2)
    stop("Neither argument for * is a functional data object fdata.")
  if (inhe.fdata1 && !inhe.fdata2) {
    fdataobj <- fdata1
    fdataobj[["data"]] <- fdata1[["data"]]*fdata2
  }
  if (!inhe.fdata1 && inhe.fdata2) {
    fdataobj <- fdata2
    fdataobj[["data"]] <- fdata1*fdata2[["data"]]
  }
  if (inhe.fdata1 && inhe.fdata2) {
    if (ncol(fdata1)!=ncol(fdata2)) stop("Error in columns dimensions")
    if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]])))  stop("Error in argvals")
    n1=nrow(fdata1)
    n2=nrow(fdata2)
    if (n1==n2){
      fdataobj <- fdata1
      fdataobj[["data"]] <- fdata1[["data"]]*fdata2[["data"]] } 
    else {
      if (n1<n2) {
        if (n2%%n1!=0) stop("Number of rows not compatible")
        fdataobj=fdata2
        fdataobj[["data"]]=fdata2[["data"]]*matrix(rep(t(fdata1$data),n2%/%n1),ncol=ncol(fdata1),byrow=TRUE)
      } else {
        if (n1%%n2!=0) stop("Number of rows not compatible")
        fdataobj=fdata1
        fdataobj[["data"]]=fdata1[["data"]]*matrix(rep(t(fdata2$data),n1%/%n2),ncol=ncol(fdata2),byrow=TRUE)
      }
    }
  }
  fdataobj
}
################################################################################

"/.fdata" <- function(fdata1,fdata2){
  inhe.fdata1 <-  inherits(fdata1, "fdata")
  inhe.fdata2 <-  inherits(fdata2, "fdata")
  if (!inhe.fdata1 && !inhe.fdata2)
    stop("Neither argument for / is a functional data object fdata.")
  if (inhe.fdata1 && !inhe.fdata2) {
    fdataobj <- fdata1
    fdataobj[["data"]] <- fdata1[["data"]]/fdata2
  }
  if (!inhe.fdata1 && inhe.fdata2) {
    fdataobj <- fdata2
    fdataobj[["data"]] <- fdata1/fdata2[["data"]]
  }
  if (inhe.fdata1 && inhe.fdata2) {
    if (ncol(fdata1)!=ncol(fdata2)) stop("Error in columns dimensions")
    if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]])))  stop("Error in argvals")
    n1=nrow(fdata1)
    n2=nrow(fdata2)
    if (n1==n2){
      fdataobj <- fdata1
      fdataobj[["data"]] <- fdata1[["data"]]/fdata2[["data"]] } 
    else {
      if (n1<n2) {
        if (n2%%n1!=0) stop("Number of rows not compatible")
        fdataobj=fdata2
        fdataobj[["data"]]=matrix(rep(t(fdata1$data),n2%/%n1),ncol=ncol(fdata1),byrow=TRUE)/fdata2[["data"]]
      } else {
        if (n1%%n2!=0) stop("Number of rows not compatible")
        fdataobj=fdata1
        fdataobj[["data"]]=fdata1[["data"]]/matrix(rep(t(fdata2$data),n1%/%n2),ncol=ncol(fdata2),byrow=TRUE)
      }
    }
  }
  fdataobj
}
################################################################################

"[.fdata"  <-  function(fdataobj, i = TRUE, j = TRUE,drop=FALSE) {
if (is.numeric(j) && j==1 && length(j)==1)
fdataobj[["data"]]  <-  matrix(fdataobj[["data"]][i,j],nrow=1)
else  fdataobj[["data"]]  <-  fdataobj[["data"]][i,j,drop=drop]
fdataobj[["argvals"]]  <-  fdataobj[["argvals"]][j]
fdataobj
}
################################################################################

"!=.fdata"  <-  function(fdata1,fdata2){
eps=1e-14
fdataequal  <-  TRUE
 if (!(all(fdata1[["data"]] == fdata2[["data"]]))) {
    res <- fdata1[["data"]] - fdata2[["data"]]
    if (!all(abs(res)<eps))        {
       fdataequal  <-  FALSE  }
#        print("No equal data matrix")
    }
 if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]]))) {
        fdataequal  <-  FALSE
        print("No equal argvals vector")
    }
 if (!(all(fdata1[["rangeval"]] == fdata2[["rangeval"]]))) {
        fdataequal  <-  FALSE
        print("No equal rangeval vector")
    }
 return(!fdataequal)
}

################################################################################

"==.fdata" <- function(fdata1,fdata2){
eps=1e-14
fdataequal <- TRUE
d1 <- dim(fdata1)
d2 <- dim(fdata2)
if (d1[1]!=d2[1]) return(FALSE)#print("Different dimensions in rows")
if (d1[2]!=d2[2]) return(FALSE)#print("Different dimensions in columns")
 if (!(all(fdata1[["data"]] == fdata2[["data"]]))) {
    res <- fdata1[["data"]] - fdata2[["data"]]
    if (!all(abs(res)<eps))        {
       fdataequal  <-  FALSE  }
#        print("No equal data matrix")
    }
 if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]]))) {
        fdataequal  <-  FALSE
        print("No equal argvals vector")
    }
 if (!(all(fdata1[["rangeval"]] == fdata2[["rangeval"]]))) {
        fdataequal  <-  FALSE
        print("No equal rangeval vector")
    }
 return(fdataequal)
}

################################################################################

"^.fdata" <- function(fdataobj,pot){
if (!inherits(fdataobj, "fdata"))  stop("No fdata object")
fdataobj[["data"]] <- fdataobj[["data"]]^pot
fdataobj
}

################################################################################

dim.fdata <- function(x) {dim(x[["data"]])}
ncol.fdata  <-  function(x){ncol(x[["data"]])}
nrow.fdata  <-  function(x){nrow(x[["data"]])}


length.fdata  <-  function(x){
  NROW(x[["data"]])
}
NROW.fdata  <-  function(x){  
  NROW(x[["data"]])
}

NCOL.fdata  <-  function(x){  
  NCOL(x[["data"]])
}



c.fdata <- function(...) {
    C=match.call()    
    fdatalist  <-  list(...)    
    n  <-  length(fdatalist)
    ii <- 2
    if (!is.null(fdatalist[[1]])) fdata1  <-  fdatalist[[1]]
    else  {
     fdata1  <-  fdatalist[[2]]
     ii <- min(3,n)
     n <- n-1
           }
    v <- rep(FALSE,len=n)

    if (n == 1)  return(fdata1)
    if (is.vector(fdata1$data))  {
         fdata1$data=matrix(fdata1$data,nrow=1)
         v[1] <- TRUE
         }
    if (is.null(fdata1)) data <- fdata1
    else data <-  fdata1$data
    dimdata  <-  dim(data)
    ndim  <-  length(dimdata)
    argvals  <-  fdata1$argvals
    rangeval  <-  fdata1$rangeval
    names  <-  fdata1$names
   if (!inherits(fdata1, "fdata"))  stop("Objects must be of class fdata")
    for (j in (ii:n)) {
        fdataj  <-  fdatalist[[j]]
        if (is.vector(fdataj$data))  {
           fdataj$data=matrix(fdataj$data,nrow=1)
           v[j] <- TRUE
           fdatalist[[j]] <- fdataj
           }
        if (!inherits(fdataj, "fdata"))
            stop("Objects must be of class fdata")                
        if (any(unlist(fdataj$argvals) != unlist(argvals)))
            stop("Objects must all have the same argvals")
        if (any(unlist(fdataj$rangeval) != unlist(rangeval)))
            stop("Objects must all have the same rangeval")
#        if (any(unlist(fdataj$names) != unlist(names)))        {
            #print("Concatenate main names")
#            names$main <- paste(names$main,"_",fdataj$names$main,sep="")
#            }
        if (length(dim(fdataj$data)) != ndim)
            stop("Objects must all have the same number of multiple functions")
    }
    if (ndim == 2) {
        for (j in ii:n) {
           fdataj  <-  fdatalist[[j]]
           dataj  <-  fdataj$data
           dd <- C[[j+1]]
           if (v[j])    rownames(dataj) <- deparse(substitute(dd))
           data  <-  rbind(data, dataj)
        }
        dd1 <- C[[2]]
       if (v[1])         rownames(data)[1] <- deparse(substitute(dd1))
    }
    concatfdata  <-  fdata(data, argvals,rangeval, names)
    return(concatfdata)
}

################################################################################

argvals <- function(fdataobj){
    if (!inherits(fdataobj, "fdata")) stop("Object must be of class fdata")
    fdataobj$argvals
}

rangeval <- function(fdataobj){
    if (!inherits(fdataobj, "fdata")) stop("Object must be of class fdata")
    fdataobj$rangeval
}



################################################################################
 "[.fdist" <- function(fdataobj, i = TRUE, j = TRUE,drop=FALSE) {
a1 <- class(fdataobj)
a2 <- attr(fdataobj,"call")
a3 <- attr(fdataobj,"par.metric")
class(fdataobj) <- "matrix"
if (is.numeric(j) && j==1 && length(j)==1) fdataobj <- matrix(fdataobj[i,j],nrow=1)
else {
#fdataobj <- fdataobj[i,j,drop=TRUE]
fdataobj <- fdataobj[i,j,drop=drop]
}
#class(fdataobj)  <-  a1
attr(fdataobj,"call") <- a2
attr(fdataobj,"par.metric") <- a3
#class(fdataobj) <- "matrix"
invisible(fdataobj)
}

################################################################################
# Deprecated
#count.na <- function(A){any(is.na(A))}
#####################
is.na.fdata <- function(x){
  if (inherits(x,"fdata")) apply(x$data,1,anyNA)
  else warning("No fda.usc class object")
}
#####################
anyNA.fdata <- function(x, recursive = FALSE){
  if (inherits(x,"fdata")) anyNA(x$data, recursive = recursive)
  else warning("No fda.usc class object")
}
count.na.fdata <- function(x){
  if (inherits(x,"fdata")) rowSums(is.na(x$data))
  else warning("No fda.usc class object")
}
##################### 
# is.square.fdata <- function( x )
# {
#   if ( is.fdata( x ) )  x <- x$data
#   if ( is.matrix( x ) ) return( nrow(x) == ncol(x) )  
#   else stop( "argument x is not a matrix or fdata class" )
# }
#####################
# is.symmetric.fdata <- function( x )
# {
#   if ( is.fdata( x ) )  x <- x$data
#   if ( !is.matrix( x ) ) {
#     stop( "argument x is not a matrix" )
#   }
#   if ( !is.numeric( x ) ) {
#     stop( "argument x is not a numeric matrix" )
#   }    
#   if ( !is.square.matrix( x ) )
# stop( "argument x is not a square numeric matrix" )
#   return( sum( x == t(x) ) == ( nrow(x) ^ 2 ) )
# }
#####################
# matrix.trace<-function (x) #tracemat<-function(x)
# {
#   if (!is.square.matrix(x)) 
#     stop("argument x is not a square matrix")
#   return(sum(diag(x)))
# }
##################### 
#DEPRECATED
traza<-function (A) 
{
  if (!is.fdata(A)) 
    sum(diag(A), na.rm = TRUE)
  else sum(diag(A$data), na.rm = TRUE)
}