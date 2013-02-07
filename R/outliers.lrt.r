outliers.lrt<-function(fdataobj,nb=200,smo=0.05,trim=0.10,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
n<-nrow(fdataobj)
    m<-ncol(fdataobj)
    if (is.null(n) && is.null(m)) stop("Error in in the data dimensions")
    out.thres.lrt<-outliers.thres.lrt(x,nb=nb,smo=smo,trim=trim,...)
    hay<-1
    outliers<-c()
    valor.estadistico<-c()
    nout<-0
    ngood<-n-nout
    curvasgood<-fdataobj
    i<-1
    maxiter<-5
    while (hay==1 && i<maxiter){
          i<-i+1
          aux<-c()
          auxmean<-func.trim.mode(curvasgood,trim=trim,...)[["data"]][1,]
          aa<-func.trimvar.mode(curvasgood,trim=trim,...)
          auxdt<-sqrt(aa[["data"]][1,])
          for (j in 1:ngood){
              aux[j]<-metric.lp(curvasgood[["data"]][j,]/auxdt,auxmean/auxdt,...)}
          maximo<-as.numeric(max(aux))
          fecha<-as.numeric(row.names(curvasgood[["data"]])[maximo==aux])
          elim<-which(maximo==aux)
          if (maximo>out.thres.lrt){
             curvasgood<-curvasgood[elim,]
             outliers<-c(outliers,fecha)
             valor.estadistico<-c(valor.estadistico,maximo)
             nout<-nout+1
             ngood<-n-nout
          }
          if (maximo<out.thres.lrt){hay<-0}
    }
if (i==maxiter) print("convergencia no lograda")
return(list("outliers"=outliers,"stat.value"=valor.estadistico,
"percentile"=out.thres.lrt))
c(outliers,valor.estadistico,out.thres.lrt)
}

