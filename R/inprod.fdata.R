################################################################################
inprod.fdata=function (fdata1,fdata2=NULL, w = 1, ...)   {
if (!is.fdata(fdata1)) stop("No fdata class")
if (is.null(fdata2)) {fdata2<-fdata1}
else  if (!is.fdata(fdata2)) stop("No fdata class")
DATA1<-fdata1[["data"]]
tt<-fdata1[["argvals"]]
DATA2<-fdata2[["data"]]
tt2<-fdata2[["argvals"]]
rtt<-fdata1[["rangeval"]]
numgr = nrow(DATA1)
numgr2 = nrow(DATA2)
n =  ncol(DATA1)
equi=TRUE
if  (sum(tt!=tt2)!=0) stop("Error: different discretization points in the input data.\n")
if (length(w) == 1) {w = rep(w, (n - 1))    }
if (length(w)!=(n-1)) {stop("DATA ERROR: The weight vector hasn't the length of the functions\n")}
 di=diff(tt);
 mdist = array(0, dim = c(numgr, numgr2))
 for (i in 1:numgr) {
    for (ii in 1:numgr2) {
         f=DATA1[i,]*DATA2[ii,]   ###
         if (equi) mdist[i,ii]=(sum((f[-n])*w)+sum((f[-1])*w))/(2*(n-1))
         else {             # tt no equiespaciados
            ai2=diff(tt)/rtt
            f2=DATA1[i,]*DATA2[ii,]
            mdist[i,ii]=sum(((f2[1:(n-1)])*ai2*w)/(2*(diff(rtt))))
}}}
mdist<-mdist*((diff(rtt)))
return(mdist)
}

