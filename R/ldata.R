#' @rdname ldata
#' @name ldata
#' @title ldata class definition and utilities  
#' 
#' @description ldata is a listt with two type of objects:
#' \itemize{
#' \item \code{df} is a data frame with the multivariate data with n rows.
#' \item \code{...}  fdata objects of class \code{fdata} with n rows.
#' }
#' @aliases ldata c.ldata  is.ldata names.ldata [.ldata plot.ldata
#' @param ldata,x object of class \code{ldata}
# @param f index of observations
# @param drop passed on to respcetive function
#' @param i index
#' @param subset subset
#' @param \dots Further arguments passed to methods.
#' @param ask logilcal    If TRUE (and the R session is interactive) the user is asked for input, before a new figure is drawn. 
#' @param color colors to interpolate; must be a valid argument to  \code{colorRampPalette}.
#' @param var.name name of continuous univariate variable used in \code{color} argument
#' @param df data frame 
#' @examples 
#' data(tecator)
#' ab0 <- tecator$absorp.fdata
#' ab1 <- fdata.deriv(ab0)
#' ab2 <- fdata.deriv(ab0,nderiv=2)
#' ldat<-ldata(tecator$y,ab1,ab2)
#' is.ldata(ldat)
#' class(ldat)
#' plot(ldat[[1]])
#' plot(ldat[[2]]) 
#' # plot(ldat)
#' # plot(ldat,var.name="Fat")
#' 

#' @export 
ldata <-function(df,...){
  C <- match.call()
  # a <- list()
  # mf <- match.call(expand.dots = FALSE)
  # m <- match("df", names(mf), 0L)
  nam <- as.character(C[3:length(C)])
  mdat <- list(...)
  clases <- sapply(mdat,class)
  if (!all(clases=="fdata")) 
    stop("ldots must be fdata objects")
  n <- sapply(mdat,NROW)
  if (!all(n==n[1])) 
    stop("Different number of rows in the fdata objects")
  if (missing(df)) 
    df <- data.frame("index"=1:n[1])
  if (!is.data.frame(df)) 
    stop("df argument must be a data frame object")
  names(mdat) <- nam
  ldat <- c(list("df"=df),mdat)
  # ldat <- list("df"=df)
  # ldat <- c(ldat,mdat)
  class(ldat) <- c("ldata","list")
  return(ldat)
}

#' @export 
c.ldata<-function (x, f, drop = FALSE, ...) 
{
  if (!is.list(x)) stop("x is not a list") 
  lenl<-length(x)
  out<-x
  lenf<-length(f)
  for (i in 1:lenl) {
    if (lenf>nrow(x[[i]])) stop("Incorrect length of f")
    out[[i]] <- x[[i]][f,]
  }
  out
}


#' @rdname ldata
#' @export 
names.ldata <- function(x){
  # x<-ldf
  if (any(class(x)=="ldata")){
    class(x)<-"list"
    nam.x <- names(x)
    nam.df  <-names(x$df)
    ind.df <- which(nam.x=="df")
    name <-list("names.df" =names(x$df),"names.fdata"=nam.x[-ind.df])
  }else stop("No ldata class object")
  return(name)
}

#' @rdname ldata
#' @export is.ldata
is.ldata=function(x){
  nc=length(x)
  nam=names(x)
  n=NROW(x[[1]])
  lo1=all(unlist(sapply(x,NROW))==n)
  if (!lo1) warning("The number of rows are different")
  clas=unlist(lapply(x,class))
  l1=which(clas=="data.frame")
#  print(l1)
 # l1=which(nam=="data.frame")
  nam<-names(clas)
  lo2=ifelse(length(l1)==0,FALSE,(length(l1)==1 & (nam[l1]=="df" | nam[l1]=="names.df")))
  lo3= ifelse(length(l1)==0,FALSE,all(clas[-l1]=="fdata"))
  if (!lo3 & !lo2) warning("The structure list(df,fdata1,fdata2,...) is not correct")
  return( lo1 & lo2 & lo3)
}

#' @rdname ldata
#' @export 
"[.ldata"=function(ldata,i){
  if (missing(i)) return(ldata)
  for (m in seq_along(ldata)){
    if (class(ldata[[m]])=="data.frame") {ldata[[m]]=ldata[[m]][i,]} else {ldata[[m]]=ldata[[m]][i]}
  }
  return(ldata)
}


#' @rdname ldata
#' @export 
subset.ldata<-function(x, subset,...){
  #if (any(class(x)!="lfdata")) stop("No list class object")
  nvar<-length(x)
  namx=names(x)
  if (missing(subset)) subset<-!logical(nrow(x[[1]]))
  if (is.numeric(subset) & max(subset)<nrow(x[[1]])) {
    subset2<-logical(nrow(x[[1]]))
    subset2[subset]<-TRUE
    subset<-subset2    
  }  
  newx<-x
  for (i in 1:nvar){
    if (is.fdata(x[[i]])) {
      newx[[i]]<-subset.fdata(x[[i]],subset,drop=FALSE,...)
    } else {
      newx[[i]]<-subset(x[[i]],subset,drop=FALSE,...)	
    }
  } 
  return(invisible(newx))
}

#' @rdname ldata
#' @export
plot.ldata <- function(x, ask=FALSE, color, var.name,...){
  if (!is.ldata(x)) stop("No ldata class object")
  #if (is.ldata) stop("No ldata class object")
  
  col.bar=FALSE
  if (!missing(var.name)){
    if (missing(color)) color=c("red","blue")
    var.color <- x$df[,var.name]
    col.bar<-TRUE
    if (!is.factor(var.color))  {
      nticks <- 5
      color=colorRampPalette(color, alpha = TRUE)(nticks)
      xfact=cut( var.color,quantile( var.color,seq(0,1,len=nticks+1)),include.lowest=TRUE)
      min.col=min(var.color)
      max.col=max(var.color)
      
    }   else {
      nticks <-nlevels(var.color)
      color=colorRampPalette(color, alpha = FALSE)(nticks)
      xfact<-var.color
      min.col=0#min(as.numeric(var.color))
      max.col=1#max(as.numeric(var.color))
      min.col=min(as.numeric(var.color))
      max.col=max(as.numeric(var.color))
    }
  }  else     {
    if (missing(color)) color=1:nrow(x$df)
    xfact<-1:nrow(x$df)
  }
  mf=5
  nvar<-length(x)-1
  if (nvar>4) ask=TRUE
  if (ask) {par(mfrow = c(1, 1))
    dev.interactive()
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }  else{    
    mf<-switch(nvar,
               "1"={c(1,1)},
               "2"={c(1,2)},
               "3"={c(1,3)},
               "4"={c(2,2)})            
    par(mfrow =mf)                    
  }
  names1<-names2<-names<-x[[1]][["names"]]
  names1$main<-"Multivariate Functional Data"
  #    tr<-paste("mode.tr",trim*100,"\u0025",sep="")         
  nam<-names(x)
  idf<-which(nam=="df")
  #plot(ts(ldf$df))
  #warning("Only fdata objects are plotted")
  nam<-nam[-idf]
  nam<-names(x)$names.fdata
  for (idat in nam) {
    data<-x[[idat]]$data
    tt<-x[[idat]]$argvals
    rtt<-x[[idat]]$rangeval
    plot(x[[idat]], col =  color[xfact],
               main  =paste0(idat,"-",x[[idat]]$names$main),...)
    
    #plot.fdata(ab,lty=1,col=colores[cfat],main="Original Trajectories")
    if (col.bar)    color.bar(color,min.col,max.col,ro=0,ticks=var.color)
    #if (col.bar)    color.bar(color,min.col,max.col,ro=0)
  }
}


################################################################################
"[.lfdata"=function(lfdata,i){
  if (missing(i)) return(lfdata)
  res=lapply(lfdata,"[",i)
  class(res)="lfdata"
  return(res)
}
################################################################################

################################################################################
is.lfdata=function(lfdata){
  n=nrow(lfdata[[1]])
  lo1=all(unlist(lapply(lfdata,nrow))==n)
  if (!lo1) warning("The number of rows are different")
  lo2=all(unlist(lapply(lfdata,is.fdata)))
  if (!lo2) warning("Not all components are fdata")
  return( lo1 & lo2 )
}
################################################################################



# @export
lfdata=function(lfdata){
n=nrow(lfdata[[1]])
lo1=all(unlist(lapply(lfdata,nrow))==n)
if (!lo1) warning("The number of rows are different")
lo2=all(unlist(lapply(lfdata,is.fdata)))
if (!lo2) warning("Not all components are fdata")
return( lo1 & lo2 )
}

############################################
# data(tecator)
# x <- tecator[[1]]
# tecator$y$cat<-factor(ifelse(tecator$y$Fat<15,2,4))
# ldf<-list("df"=tecator$y,"x"=x,"x.d1"=fdata.deriv(x))
# plot.mfdata(ldf,col=list("x"=2,"x.d1"=4))
# plot.mfdata(ldf,col=2:3)

# plot.lfdata<-function(lfdata,ask=FALSE,color,...){
#   mf=5
#   nvar<-length(lfdata)
#   if (nvar>4) ask=TRUE
#   if (ask) {par(mfrow = c(1, 1))
#     dev.interactive()
#     oask <- devAskNewPage(TRUE)
#     on.exit(devAskNewPage(oask))}
#   else{    mf<-switch(nvar,
#                       "1"={c(1,1)},
#                       "2"={c(1,2)},
#                       "3"={c(1,3)},
#                       "4"={c(2,2)})            
#   par(mfrow =mf)                    }
#   names1<-names2<-names<-lfdata[[1]][["names"]]
#   names1$main<-"Multivariate Functional Data"
#   #    tr<-paste("mode.tr",trim*100,"\u0025",sep="")         
#   nam<-names(lfdata)
#   
#   if (is.null(nam)) nam<-1:nvar
#   for (idat in 1:nvar) {
#     data<-lfdata[[idat]]$data
#     tt<-lfdata[[idat]]$argvals
#     rtt<-lfdata[[idat]]$rangeval
#     if (missing(color)) color2<-1
#     else {
#       if (is.list(color)) color2<-color[[idat]]
#       else color2<-color
#     }
#     plot(lfdata[[idat]], col =  color2,lty=1, main =nam[idat],...)
#   }
#   
# }    

