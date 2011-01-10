plot.fdata<-function(x,type,main,xlab,ylab,...) {
if (!is.fdata(x))  stop("Object is not fdata class")
if (missing(type)) type="l"
if (missing(main)) main=x[["names"]][["main"]]
if (missing(xlab)) xlab=x[["names"]][["xlab"]]
if (missing(ylab)) ylab=x[["names"]][["ylab"]]
if (is.vector(x[["data"]])) matplot(x[["argvals"]],(x[["data"]]),type=type,main=main,ylab=ylab,xlab=xlab,...)
else matplot(x[["argvals"]],t(x[["data"]]),type=type,main=main,ylab=ylab,xlab=xlab,...)
}



#lines.fdata<-function(x,...) {
#if (!is.fdata(x))  stop("Object is not fdata class")
#lines(x[["argvals"]],t(x[["data"]]),...)}


lines.fdata=function(x,...){plot(x,add=TRUE,...)}


title.fdata<-function(x,main=NULL,xlab=NULL,ylab=NULL,rownames=NULL) {
if (!is.fdata(x))  stop("Object is not fdata class")
if (!is.null(rownames)) rownames(x[["data"]])<-rownames
if (!is.null(main)) x[["names"]][["main"]]<-main
if (!is.null(xlab)) x[["names"]][["xlab"]]<-xlab
if (!is.null(ylab)) x[["names"]][["ylab"]]<-ylab
x
}
