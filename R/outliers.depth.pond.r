outliers.depth.pond<-function(fdataobj,nb=200,smo=0.05,dfunc=depth.mode,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
n<-nrow(fdataobj)
m<-ncol(fdataobj)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
cutoff<-median(quantile.outliers.pond(x,dfunc=dfunc,nb=nb,smo=smo,...))
    hay<-1
    outliers<-c();dep.out<-c()
    curvasgood<-fdataobj
    row.names(curvasgood[["data"]])=1:n
    while (hay==1){
          d=dfunc(curvasgood,...)$dep
          if (is.null(outliers)){dtotal<-d}
          fecha<-as.numeric(row.names(curvasgood[["data"]])[d<cutoff])
          elim<-which(d<cutoff)
          if (length(elim)>0){
             dep.out<-c(dep.out,d[d<cutoff])
             curvasgood<-curvasgood[-elim,]
             outliers<-c(outliers,fecha)
          }
        if (length(elim)==0 || length(outliers)>n/5){hay<-0}
    }
return(list("outliers"=outliers,"quantile"=cutoff,"dep"=dtotal,"dep.out"=dep.out))
c(outliers,cutoff,dtotal,dep.out)
}



