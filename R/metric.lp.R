metric.lp=function (fdata1, fdata2 = NULL, lp = 2, w = 1, ...)
{
p<-lp
C1<-match.call()
same<-FALSE
if (is.fdata(fdata1)) {
   fdat<-TRUE
   tt<-tt1<-fdata1[["argvals"]]
   rtt<-fdata1[["rangeval"]]
   nas1<-apply(fdata1$data,1,count.na)
   if (any(nas1))  stop("fdata1 contain ",sum(nas1)," curves with some NA value \n")
   if (is.null(fdata2)) {fdata2<-fdata1;same<-TRUE}
   else  if (!is.fdata(fdata2))  {

         fdata2<-fdata(fdata2,tt1,rtt,fdata1$names)
         }
   nas2<-apply(fdata2$data,1,count.na)
   if (any(nas2)) {stop("fdata2 contain ",sum(nas2)," curves with some NA value \n")}
   DATA1<-fdata1[["data"]]
   DATA2<-fdata2[["data"]]
   tt2<-fdata2[["argvals"]]
   if (!same) {
   if  (sum(tt1!=tt2)!=0) {stop("Error: different discretization points in the input data.\n")  }
   }}
else{
    fdat<-FALSE
     DATA1<-fdata1
     if (is.null(fdata2)) {fdata2<-fdata1;same<-TRUE}
     DATA2<-fdata2
}
numgr = nrow(DATA1)
numgr2 = nrow(DATA2)
np<-ncol(DATA1)
if ((length(w) != np) & (length(w) != 1)) {
   stop("DATA ERROR: The weight vector hasn't the length of the functions\n")
}
testfordim <- sum(dim(DATA1) == dim(DATA2)) == 2
twodatasets <- TRUE
if (testfordim)     twodatasets <- sum(DATA1-DATA2,na.rm=TRUE)==0
if (fdat) {
 dtt<-diff(tt)
 eps<-as.double(.Machine[[1]]*100)
 inf<-dtt-eps;sup<-dtt+eps
 if (all(dtt>inf) & all(dtt<sup)) {equi=TRUE}
 else equi=FALSE
 mdist = array(0, dim = c(numgr, numgr2))
 predi <- TRUE
 if (testfordim) {
   if (twodatasets) {
        predi <- FALSE
        for (i in 1:numgr) {
           ii = i + 1
           for (ii in i:numgr2) {
             f=w*abs(DATA1[i,]-DATA2[ii,])^p
           mdist[i,ii]=(int.simpson2(tt,f,equi))^(1/p)
                   }      }
    mdist = t(mdist) + mdist
}
#mdist<-mdist/(sqrt(rtt[2]-rtt[1]))
}
if (predi) {
   for (i in 1:numgr) {
     for (ii in 1:numgr2) {
             f=w*abs(DATA1[i,]-DATA2[ii,])^p
                  mdist[i,ii]=(int.simpson2(tt,f,equi))^(1/p)
             }}
#mdist<-mdist/(sqrt(diff(rtt)))
}
}
else {
   mdist = array(0, dim = c(numgr, numgr2))
   for (i in 1:numgr) {
     for (ii in 1:numgr2) {
            f=w*(abs(DATA1[i,]-DATA2[ii,])^p)
            mdist[i,ii]=(sum(f))^(1/p)
             }}
}
attr(mdist,"call")<-"metric.lp"
attr(mdist,"par.metric")<-list("lp"=p,"w"=w)
return(mdist)
}

