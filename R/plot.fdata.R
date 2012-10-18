plot.fdata<-function(x,type,main,xlab,ylab,...) {
if (any(class(x)=="fdata2d"))  {
#stop("Object is not fdata2d class")
if (missing(type)) type="persp"
#if (missing(main)) main=x[["names"]][["main"]]
#if (missing(xlab)) xlab=x[["names"]][["xlab"]]
#if (missing(ylab)) ylab=x[["names"]][["ylab"]]

len.dm<-length(dim(x$data))
if (len.dm==2) { 
switch (type,
"persp"={                          
par(bg = "white")
xx <- x[["argvals"]][[1]]
y <- x[["argvals"]][[2]]
z <- x[["data"]]
#nrz <- nrow(z);
#ncz <- ncol(z)
#jet.colors <- colorRampPalette( c("yellow", "red") ) 
#nbcol <- length(xx)
#color <- jet.colors(nbcol)
#zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
#facetcol <- cut(zfacet, nbcol)
 persp(x=xx,y=y,z=z,xlim=x[["rangeval"]][[1]],ylim=x[["rangeval"]][[2]],
 main =  x$names$main,xlab = x$names$xlab[1], ylab =  x$names$xlab[2],...)
#, col=color[facetcol],...)
# par(op)
},
"filled.contour"={
filled.contour(x=x[["argvals"]][[1]],y=x[["argvals"]][[2]],z=x[["data"]],
xlim=x[["rangeval"]][[1]],ylim=x[["rangeval"]][[2]],
plot.title=title(main = x$names$main,
xlab =x$names$xlab[1], ylab =  x$names$xlab[2]),...)},

"contour"={contour(x=x[["argvals"]][[1]],y=x[["argvals"]][[2]],z=x[["data"]],
xlim=x[["rangeval"]][[1]],ylim=x[["rangeval"]][[2]],
plot.title=title(main =  x$names$main,
    xlab = x$names$xlab[1], ylab =  x$names$xlab[2]),...)},#labels repetidos
    
"image"={image(x = x[["argvals"]][[1]],y = x[["argvals"]][[2]],z= x[["data"]],
xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],main =  x$names$main[1],
    xlab = x$names$xlab[1], ylab =  x$names$xlab[2],...)},
    
"contourplot"={contourplot(data=x[["data"]],
xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],...)}    
)
}
else {
  if (len.dm==3) { contourplot(data=x[["data"]],
xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],...)}
    }
if (len.dm>3) stop("Not implemented plot for arrays of more than 3 dimension yet")

}
else {
if (!is.fdata(x))  stop("Object is not fdata class")
if (missing(type)) type="l"
if (missing(main)) main=x[["names"]][["main"]]
if (missing(xlab)) xlab=x[["names"]][["xlab"]]
if (missing(ylab)) ylab=x[["names"]][["ylab"]]
if (is.vector(x[["data"]])) matplot(x[["argvals"]],(x[["data"]]),type=type,main=main,ylab=ylab,xlab=xlab,...)
else matplot(x[["argvals"]],t(x[["data"]]),type=type,main=main,ylab=ylab,xlab=xlab,...)
}
}

lines.fdata=function(x,...){plot(x,add=TRUE,...)}


title.fdata<-function(x,main=NULL,xlab=NULL,ylab=NULL,rownames=NULL) {
if (!is.fdata(x))  stop("Object is not fdata class")
if (!is.null(rownames)) rownames(x[["data"]])<-rownames
if (!is.null(main)) x[["names"]][["main"]]<-main
if (!is.null(xlab)) x[["names"]][["xlab"]]<-xlab
if (!is.null(ylab)) x[["names"]][["ylab"]]<-ylab
x
}
