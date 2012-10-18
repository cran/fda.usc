outliers.depth.pond<-function(fdataobj,nb=200,smo=0.05,dfunc=depth.mode,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
 nas1<-apply(fdataobj$data,1,count.na)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
 x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
n<-nrow(fdataobj)
m<-ncol(fdataobj)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
scale=FALSE
cutoff<-median(quantile.outliers.pond(x,dfunc=dfunc,nb=nb,smo=smo,...))
    hay<-1
    outliers<-c();dep.out<-c()
    curvasgood<-fdataobj
    row.names(curvasgood[["data"]])=1:n
    dd<-dfunc(curvasgood,...)     
    modal<-FALSE
    if (!is.null(dd$dist)) {
      modal=TRUE   
      dd<-dfunc(curvasgood,...)
          }
    d<-dd$dep 
    
    while (hay==1){
#          d=dfunc(curvasgood,...)$dep
          if (is.null(outliers)){dtotal<-d}
          fecha<-as.numeric(row.names(curvasgood[["data"]])[d<cutoff])
          elim<-which(d<cutoff)
          if (length(elim)>0){
             dep.out<-c(dep.out,d[d<cutoff])
             curvasgood<-curvasgood[-elim,]
             outliers<-c(outliers,fecha)
          }
        if (length(elim)==0 || length(outliers)>n/5){hay<-0}
                else {
            if (modal) {
             mdist<-dd$dist[-elim,-elim]
            class(mdist)<-c("matrix","fdist")        
            dd<-dfunc(curvasgood,metric=mdist,...)
            }
          else dd<-dfunc(curvasgood,...)
          d<-dd$dep 
          }
    }
outliers<-rownames(fdataobj[["data"]])[outliers]    
return(list("outliers"=outliers,"quantile"=cutoff,"dep"=dtotal,"dep.out"=dep.out))
c(outliers,cutoff,dtotal,dep.out)
}

 