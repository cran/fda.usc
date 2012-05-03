fdata=function(mdata,argvals=NULL,rangeval=NULL,names=NULL){
out=list("data"=NULL)
if (length(class(mdata))>1) class(mdata)<-class(mdata)[1]
out<-switch(class(mdata),
matrix={
        out[["data"]]=mdata
        out},
data.frame={
            out[["data"]]=as.matrix(mdata)
            out},
fdata=stop("no fue posible la conversion"),
numeric={
         out[["data"]]=matrix(mdata,nrow=1)
         out},
integer={
         out[["data"]]=matrix(mdata,nrow=1)
         out},
fd={
   r = mdata$basis[[3]]
   if (is.null(argvals))
     argvals= seq(r[1], r[2], len =mdata$basis$nbasis)
    nb<- length(argvals)
    tt = argvals
#   tt = seq(r[1], r[2], len = length(mdata$fdnames$time))

   out[["data"]] = t(eval.fd(tt, mdata))
   if (!is.null(mdata$fdnames$reps)) rownames(out[["data"]]) = mdata$fdnames$reps
   else      rownames(out[["data"]]) =1:nrow( out[["data"]])
   if (!is.null(mdata$fdnames$time)) {
      colnames(out[["data"]]) = 1:ncol( out[["data"]])
      }
   else      {   colnames(out[["data"]]) =1:ncol( out[["data"]])   }
   out
   },
fds={
 out[["data"]] = mdata$y
 if (is.null(mdata$time))       out[["argvals"]] =1:ncol(out[["data"]])
 else out[["argvals"]]= seq(mdata$time[1], mdata$time[length(mdata$time)],
  len = length(mdata$time))
  out},
fts={
 out[["data"]] = mdata$y
 if (is.null(mdata$time))       out[["argvals"]]<-1:ncol(out[["data"]])
 else out[["argvals"]]<- seq(mdata$time[1], mdata$time[length(mdata$time)],
  len = length(mdata$time))
  out},
sfts={
 out[["data"]] = mdata$y
 if (is.null(mdata$time))       out[["argvals"]]=1:ncol(out[["data"]])
 else out[["argvals"]]= seq(mdata$time[1], mdata$time[length(mdata$time)],
  len = length(mdata$time))
  out}
)
nc<-nrow(out[["data"]])
np<-ncol(out[["data"]])
if (is.null(argvals)) {
 if (is.null(colnames(out[["data"]]))) {out[["argvals"]]=1:ncol(out[["data"]])}
   else {
#   if (!any(is.na(as.numeric(colnames(out[["data"]]))))) {
#    out[["argvals"]]=as.numeric(colnames(out[["data"]]))   }
#    else    out[["argvals"]]=1:ncol(out[["data"]])}   }
	out[["argvals"]]=1:ncol(out[["data"]])}   }
else     out[["argvals"]]=argvals
lentt=length(out[["argvals"]])
if (is.null(rangeval)) rangeval=range(out[["argvals"]])
out[["rangeval"]]<-rangeval
if ((np!=lentt) && (nc==lentt)) {
         out[["data"]]=matrix(out[["data"]],ncol=nc)
         nc<-1
         print("Warning: The transposed data is returned")      }
else    out[["data"]]=out[["data"]]
if (is.null(dimnames(mdata))) {
#rownames(out[["data"]])<-1:nc
colnames(out[["data"]])=round(out[["argvals"]],4)
}
out[["names"]]<-list("main"="fdataobj","xlab"="t","ylab"="X(t)")
if (!is.null(names$main)) out$names$main<-names$main
if (!is.null(names$xlab)) out$names$xlab<-names$xlab
if (!is.null(names$ylab)) out$names$ylab<-names$ylab
class(out)="fdata"
return(out)
}


