metric.lp=function (fdata1, fdata2 = NULL, p = 2, w = 1, ...)
{
same<-FALSE
if (!is.fdata(fdata1)) {fdata1<-fdata(fdata1)}
DATA1<-fdata1[["data"]]
tt1<-fdata1[["argvals"]]
if (is.null(fdata2)) {fdata2<-fdata1;same<-TRUE}
else  if (!is.fdata(fdata2))  {fdata2<-fdata(fdata2) }
DATA2<-fdata2[["data"]]
tt2<-fdata2[["argvals"]]
rtt<-fdata1[["rangeval"]]
tt<-tt1
testfordim <- sum(dim(DATA1) == dim(DATA2)) == 2
twodatasets <- TRUE
if (testfordim)     twodatasets <- sum(DATA1-DATA2,na.rm=TRUE)==0
numgr = nrow(DATA1)
numgr2 = nrow(DATA2)
n =  ncol(DATA1)
equi=TRUE
if (!same) {
if  (sum(tt1!=tt2)!=0) {
 stop("Error: different discretization points in the input data.\n")  }
}
if (length(w) == 1) {w = rep(w, (n - 1))    }
if (length(w) != (n - 1)) {
      stop("DATA ERROR: The weight vector hasn't the length of the functions\n")
    }
           di=diff(tt1)
           di2=diff(tt2)
#           if (sum(di2==di)!=(n-1)) equi=FALSE
mdist = array(0, dim = c(numgr, numgr2))
predi <- TRUE
if (testfordim) {
   if (twodatasets) {
#print("no prediction")
        predi <- FALSE
        for (i in 1:numgr) {
           ii = i + 1
           for (ii in i:numgr2) {
             f=(DATA1[i,]-DATA2[ii,])
             if (equi) {
                   mdist[i,ii]=((sum((abs(f[-n])^p)*w)
                       +sum((abs(f[-1])^p)*w))/(2*(n-1)))^(1/p)
             if (sum(w)==(n-1))  mdist[i,ii]=(sum(abs(f)^p))^(1/p)/2 #np?
                       }
              else {
                 ai2=diff(tt)/rtt
                 f2=abs(DATA1[i,]-DATA2[ii,])
                 mdist[i,ii]=sum((((f2[1:(n-1)]^p)*ai2*w)/(2*(diff(rtt))))^(1/p))
        }        }      }
    mdist = t(mdist) + mdist
}
mdist<-mdist#/(sqrt(rtt[2]-rtt[1]))
}
if (predi) {
#print("prediction")
   for (i in 1:numgr) {
     for (ii in 1:numgr2) {
         f=(DATA1[i,]-DATA2[ii,])
         if (equi) mdist[i,ii]=((sum((abs(f[-n])^p)*w)+sum((abs(f[-1])^p)*w))/(2*(n-1)))^(1/p)
         else {
            ai2=diff(tt)/rtt
            f2=abs(DATA1[i,]-DATA2[ii,])
            mdist[i,ii]=sum((((f2[1:(n-1)]^p)*ai2*w)/(2*(diff(rtt))))^(1/p))
}}}
mdist<-mdist*(sqrt(diff(rtt))) ##produce NAs
}
#print("end metric.lp")
return(mdist)
}


