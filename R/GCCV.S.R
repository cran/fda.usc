GCCV.S=function(y,S,criteria="GCCV1",W=NULL,trim=0,draw=FALSE,metric=metric.lp,...){
    isfdata<-is.fdata(y)
#    tab=list("GCV","AIC","FPE","Shibata","Rice")
    tab=list("GCCV1","GCCV2","GCCV3","GCV")
    type.i=pmatch(criteria,tab)
    n=ncol(S);l=1:n

    if (isfdata) #to smooth function (min.np,min.basis)
    { 
         nn<-ncol(y)
         if (is.null(W)) W<-diag(nn)
         y2=t(y$data)
         y.est=t(S%*%y2)
         y.est<-fdata(y.est,y$argvals, y$rangeval, y$names)
         e <- y - y.est
#         e$data<-sqrt(W)%*%(e$data)   
         ee <- drop(norm.fdata(e,metric=metric,...)[,1]^2)
         if (trim>0) {
            e.trunc=quantile(ee,probs=(1-trim),na.rm=TRUE,type=4)
            ind<-ee<=e.trunc
            if (draw)  plot(y,col=(2-ind))
            l<-which(abs(ee)<=e.trunc)
            res = mean(ee[ind],na.rm=TRUE)
            }
        else  res = mean(ee, na.rm = TRUE)
        ee<-e$data
   }
   else   #to regression function
   {
        if (is.null(W)) W<-diag(n)         
        if (is.matrix(y)&&(ncol(y)==1) ){y2<-y;draw<-FALSE}
        else if (is.vector(y)){y2<-y;draw<-FALSE}
        else stop("y is not a fdata,  vector or matrix")
     y.est=S%*%y2
     e=y2-y.est
     if (trim>0) {
             ee = t(e)
             e.trunc=quantile(abs(ee),probs=(1-trim),na.rm=TRUE,type=4)
             l<-which(abs(ee)<=e.trunc)
             e<-e[l]
             }
             ee<-e
    }
    d<-diag(S)[l] 
    vv<-switch(type.i,
                   "1"={
                   Sigma<-solve(W)
                   SC<-S%*%Sigma
                   df<-traza(2*SC-SC%*%t(S))
                   mean(ee^2)/(1-df/n)^2                                                    
                   },
                   "2"={
                   Sigma<-solve(W)
                   df<-traza(S%*%Sigma)
                   mean(ee^2)/(1-df/n)^2
                   },                                                     
                   "3"={
                   Sigma<-solve(W)
                   df<-traza(S%*%Sigma%*%t(S))
                   mean(ee^2)/(1-df/n)^2
                   }
                   ,"4"={
#                   Sigma<-solve(W)
                   df<-traza(S)
                   mean(ee^2)/(1-df/n)^2
                   })
 attr(vv, "df") <- df    
 return(vv)
 }
  

################################################################################
################################################################################
GCCV.S.version.simple=function(y,S,criteria="GCV",W=diag(ncol(S)),Sigma=W,trim = 0, draw = FALSE, metric = metric.lp, 
                ...) {
  # print("entra GCV.GLSS")
  #    tab=list("GCV","AIC","FPE","Shibata","Rice","GCCV1","GCCV2","GCCV3")
  #    type.i=pmatch(criteria,tab)
  # print(type.i)    
  n=ncol(S);l=1:n
  y.est=(S%*%y)
  e <- y - y.est
  # print("entra 2")  
  vv<-switch(criteria,                  
             "GCV"=n*traza(t(e) %*% W %*% e)/(n-traza(S))^2,
             "GCCV1"={
               SC<-S%*%Sigma
               mean(e^2)/(1-traza(2*SC-SC%*%t(S))/n)^2                                                    
             },
             "GCCV2"=mean(e^2)/(1-traza(S%*%Sigma)/n)^2,                                                     
             "GCCV3"=mean(e^2)/(1-traza(S%*%Sigma%*%t(S))/n)^2                                                     
  )
  return(vv)
}
################################################################################
