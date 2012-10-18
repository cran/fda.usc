outliers.depth.trim<-function(fdataobj,nb=200,smo=0.05,trim=0.01,
dfunc=depth.mode,...){
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  nas1<-apply(fdataobj$data,1,count.na)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")         
 x<-fdataobj[["data"]]
 tt<-fdataobj[["argvals"]]
 rtt<-fdataobj[["rangeval"]]
 n<-nrow(fdataobj)
 m<-ncol(fdataobj)
# print(dfunc=depth.mode)
 if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
    cutoff<-median(quantile.outliers.trim(fdataobj,dfunc=dfunc,nb=nb,smo=smo,
    trim=trim,scale=FALSE,...))
    hay<-1
    outliers<-c()
    dep.out<-c()
    curvasgood<-fdataobj
    dd<-dfunc(curvasgood,trim=trim,...)     
    modal<-FALSE
    if (!is.null(dd$dist)) {
      modal=TRUE   
      dd<-dfunc(curvasgood,trim=trim,...)
          }
    d<-dd$dep          
    rownames(curvasgood[["data"]])=1:n
    while (hay==1){            
          if (is.null(outliers)){dtotal<-d}
          fecha<-as.numeric(rownames(curvasgood[["data"]])[d<cutoff])            
          elim<-which(d<cutoff)
          if (length(elim)>0){
             dep.out<-c(dep.out,d[d<cutoff])
             curvasgood<-curvasgood[-elim,]
             outliers<-c(outliers,elim)
      
          }    
        if (length(elim)==0 || length(outliers)>n/5){hay<-0}
        else {
            if (modal) {
             mdist<-dd$dist[-elim,-elim]
            class(mdist)<-c("matrix","fdist")        
            dd<-dfunc(curvasgood,trim=trim,metric=mdist,scale=FALSE,...)
            }
          else dd<-dfunc(curvasgood,trim=trim,...)
          d<-dd$dep 
          }
    }
outliers<-rownames(fdataobj[["data"]])[outliers]    
return(list("outliers"=outliers,"quantile"=cutoff,"Dep"=dtotal,"dep.out"=dep.out))
c(outliers,cutoff,dtotal,dep.out)
}
