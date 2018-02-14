
predict.fregre.gls<-function(object, newx = NULL, type = "response",
    se.fit= FALSE, scale = NULL,df , interval = "none", 
    level = 0.95, weights = 1, pred.var, n.ahead = 1, ...){
# print("entra p predcit f.gls")    
 if (is.null(object)) stop("No fregre.lm object entered")
 if (is.null(newx)) {
    yp=predict(object,type=type,se.fit=se.fit,interval=interval,level=level,weights=weights,pred.var=pred.var,df=df,scale=scale,...)    
    print("No newx entered")
    return(yp)
    } 
 else {
 data=newx
 basis.x=object$basis.x
 basis.b=object$basis.b
 formula=object$formula.ini
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)   
##########
 vtab<-rownames(attr(tf,"factors"))
 vnf=intersect(terms,names(data$df))
 vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
 vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 beta.l=list()
 kterms=1

# print(attributes(tf))
# print(tf)
#print("entra p predcit f.gls33")    
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
if (attr(tf,"intercept")==0) {
     print("No intecept")
     pf<- paste(pf,-1,sep="")
     }
if (length(vnf)>0) {
 first=FALSE
# print(paste("no functional variable",vnf[i]))
   for ( i in 1:length(vnf)){
#     print(paste("Non functional covariate:",vnf[i]))
     if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
#     else pf <- paste(pf, terms[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
if   (attr(tf,"intercept")==0) {
#     print("No intecept")
     pf<- paste(pf,-1,sep="")
     }
   mf<-as.data.frame(model.matrix(formula(pf),data$df))
   vnf2<-names(mf)[-1]   
    for ( i in 1:length(vnf2))  pf<-paste(pf, "+", vnf2[i], sep = "")
    XX <- mf 
}
else {
    pf2<-paste(pf, "1", sep = "")
    XX <- data.frame(model.matrix(formula(pf2),data$df))
#    print(XX)
#    names(XX)[3]<-"(Intercept)"
    first=FALSE
}
# print(names(object$XX))
#print(names(XX))
#print(vnf)
#print(dim(XX))

#yp<-object$coefficients[1]*rep(1,len=nrow(newx[[vfunc[i]]]))
if (length(vnf)>0) {
#  print(object$coefficients[names(XX)])
#  print(dim(XX))
  spm<-matrix(object$coefficients[names(XX)],ncol=1)
#  print(spm )
 yp<-as.matrix(XX)%*%spm
 }
else yp<-object$coefficients[1]*rep(1,len=nrow(newx[[vfunc[1]]])) # yp es el intercept
#print(yp)

if (length(vfunc)>0)  {
#   yp2<-a1 <- object$coefficients[1] * rep(1, len = nrow(data[[vfunc[1]]]))
   for (i in 1:length(vfunc)) {
#print(i)
# print(object$basis.x[[vfunc[i]]])
# print(object$basis.x[[vfunc[i]]]$type)
#print(vfunc)
   if(class(data[[vfunc[i]]])[1]=="fdata")  {
     fdataobj<-data[[vfunc[i]]]
      x.fd<-fdataobj[["data"]]
      tt<-fdataobj[["argvals"]]
      rtt<-fdataobj[["rangeval"]]
# print(object$basis.x[[vfunc[i]]]$type)
      if (!object$basis.x[[vfunc[i]]]$type=="pc"&!object$basis.x[[vfunc[i]]]$type=="pls") { 
#print("fda basis")                 # si es pls hay que quitarle la norma!!!
 	      	x.fd = Data2fd(argvals = tt, y = t(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]]$data),
                      basisobj = basis.x[[vfunc[i]]],fdnames=rownames(x.fd))
	    	  r=x.fd[[2]][[3]]
          J<-object$JJ[[vfunc[i]]]
          Z = t(x.fd$coefs) %*% J
#        colnames(Z) = paste(vfunc[i], ".",colnames(J), sep = "")
          colnames(Z) = colnames(J)
#print("colnames(Z)")          
      }
      else {
# print("pc o pls basis")        
          name.coef<-paste(vfunc[i], ".",rownames(object$basis.x[[vfunc[i]]]$basis$data),sep ="")
          newXcen<-fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]]                  
                      if (object$basis.x[[vfunc[i]]]$type == "pls") {
                       if (object$basis.x[[vfunc[i]]]$norm)  {
                         sd.X <- sqrt(apply(object$basis.x[[vfunc[i]]]$fdataobj$data, 2, var))
                         newXcen$data<- newXcen$data/(rep(1, nrow(newXcen)) %*% t(sd.X))
                        }
                      } 
                    Z<- inprod.fdata(newXcen,object$vs.list[[vfunc[i]]]) 
                    
#          Z<- inprod.fdata(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]],object$vs.list[[vfunc[i]]])
          colnames(Z)<-name.coef
#         object$beta.l[[vfunc[i]]]$data <- matrix(object$beta.l[[vfunc[i]]]$data,nrow = 1)
#         b1 <- inprod.fdata(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]],object$beta.l[[vfunc[i]]])
#         yp2<-yp2+b1
      }
       if (first) {    XX=Z;              first=FALSE         }
       else XX = cbind(XX, Z)
      }
      else {
          if(class(data[[vfunc[i]]])[1]=="fd")  {
             if (class(object$basis.x[[vfunc[i]]])!="pca.fd") {
             x.fd<-fdataobj<-data[[vfunc[i]]]
 	    	     r=x.fd[[2]][[3]]
             J<-object$JJ[[vfunc[i]]]
             x.fd$coefs<-x.fd$coefs-object$mean[[vfunc[i]]]$coefs[,1]
             Z = t(x.fd$coefs) %*% J
             colnames(Z) = colnames(J)
             }
             else {
                       name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(object$basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
                       data[[vfunc[i]]]$coefs<- sweep(data[[vfunc[i]]]$coefs,1,(object$basis.x[[vfunc[i]]]$meanfd$coefs),FUN="-")
                       fd.cen<-data[[vfunc[i]]]
                       #          fd.cen<-data[[vfunc[i]]]-object$basis.x[[vfunc[i]]]$meanfd # de las CP basi
                       Z<- inprod(fd.cen,object$basis.x[[vfunc[i]]]$harmonics)
                       colnames(Z)<-name.coef[[vfunc[i]]]
                   }
              if (first) {    XX=Z;              first=FALSE         }
             else XX = cbind(XX, Z)
           }
          else stop("Please, enter functional covariate")
       }  }
       }
#print(object$rn)
 n.ahead<-nn<-nrow(XX) 
 n<-length(object$residuals) 
 if (!is.data.frame(XX)) XX=data.frame(XX)
#print(object$rn)
  if (!object$rn)  {
#    print(" predict")
#    print("aa2")
    if (class(object)[1]=="gls"){  
# print(names(XX))    
# print(names(object$XX))    
     if (is.null(object$correlation)) return(predict(object=object,newdata=XX)  )    #call to nlme:::predict.gls            
# realmente hay que buscar el form Y PONER EL INDICADOR
#    print(" buscar en correlation el AR y hacer el predict")
    ype<-numeric(nn)
    rho<-coef(object$modelStruct[[1]],FALSE)    
#    XX=data.frame(XX,newx[["df"]])
    XX=data.frame(XX,newx$df)


# print(se.fit)
# print("0")
# oo<-lm(object$formula,data=object$XX)
# # ypx<-predict.lm(object=oo,newdata=object$XX[1:10,],se.fit=se.fit)    
# print(ypx)
# print("1")
# ypx<-nlme::predict.gls(object=object,newdata=object$XX[1:10,])    
# # print(ypx)
#print("2")
ypx<-predict.gls(object=object,newdata=XX)     #call to internal function predict.gls an not to nlme:::predict.gls 
#print("3")
# print(ypx)
# print(ypx)    
#    ypx2<-predict.lm(object=object,newdata=XX)    
# print(ypx)    

#class(res3$correlation)[1]    ==corARMA o corAR1
    yp<-switch(class(object$correlation)[1],
     "corAR1"={  
       #      print("corAR111")
#    if (names(rho)!="range"){        
    if (!is.null(object$correlation))       groups <- nlme::getGroupsFormula(object$correlation)
    else groups <- NULL
    cls1<-nlme::getCovariate(object$correlation,data=data$df)
    if (is.null(groups)){ 
      ype[1] <- rho *object$residuals[n]      
      if (nn>1){
      for (i in 2:nn)   ype[i] <- rho * ype[i-1] 
     }
       yp<-ypx+ype
     }
     else{
#      warning("Conditional predictions for each group are not yet done.")
    gr<-  nlme::getGroups(object$correlation,data=data$df)   
    ind<- nlme::getCovariate(object$correlation,data=data$df)
    ind.old<- nlme::getCovariate(object$correlation,data=object$data$df)
    gr.old<-  nlme::getGroups(object$correlation,data=object$data$df)    
#    (getCovariateFormula(res2$correlation)[[2]][[3]])

#    print(gr);    print("gr")
#    print(ind)    ;    print("ind")
#    print(gr.old);    print("gr.old")
#    print(ind.old)    ;    print("ind.old")    
#      gr<-(data[["df"]][,groups[[2]]])
      lengr<-nlevels(gr)
      name.group<-levels(gr)     
      for (i in 1:length(cls1)){
#cat(" entra grupo ",i)      
            ii<-gr.old==name.group[i]
            ii2<-gr==name.group[i]            
            ng<-sum(ii)#o name.group[i]           
            ype2<-numeric(ng)        
            order.ind.old<-order(ind.old[[i]])
            order.ind<-order(ind[[i]])
            res.gr<-object$residuals[ii][order.ind.old]     
            ype2[1] <- rho *res.gr[ng]             
            if (ng>1){
              for (j in 2:ng)   ype2[j] <- rho * ype2[j-1] 
              ype2<-ype2[order.ind]
            }
#print(length(ype2))            
#print(length(ype))        
#print(length(ii2))                
#print(sum(ii2))            
     ype[ii2]<-ype2
     }
#alternativa al bucle
# print("AR1")
# print(ype)
ab<-arima(object$residuals, order = c(1,0,0), fixed=coef(object$modelStruct,unconstrained=FALSE)
          ,include.mean=FALSE)
ab<-predict(ab,n.ahead=n.ahead)
print(ab)

yp<-ypx+ab$pred
# print(yp)
# print("yp")
#     yp<-ypx+ype
     }
     },
    "corARMA"={  
      #print("corARMA")
#print(rho)
    if (!is.null(object$correlation))       groups <- nlme::getGroupsFormula(object$correlation)
    else groups <- NULL
    cls1<-nlme::getCovariate(object$correlation,data=data$df)  
    p<-attributes(object$modelStruct[[1]])$p
    q<-attributes(object$modelStruct[[1]])$q 
    gr<-  nlme::getGroups(object$correlation,data=data$df)   
    ind<- nlme::getCovariate(object$correlation,data=data$df)
    ind.old<- nlme::getCovariate(object$correlation,data=object$data$df)
    gr.old<-  nlme::getGroups(object$correlation,data=object$data$df)    
#    if (is.list(ind)) ya esta bien,,sino tengo que coger el ultimo residuo del grupo y hacer el invento
#    print(gr);    print("gr")
#    print(ind)    ;    print("ind")
#    print(gr.old);    print("gr.old")
#    print(ind.old)    ;    print("ind.old")     
    rho1<-rho2<-NULL
# print(p)    
# print(q)
    if (p>0) { 
      rho1<-rho[1:p]
      if (q>0) rho2<-rho[(p+1):(p+q)]
      }
    else rho2<-rho      
# print(rho1)    
    maxpq<-max(p,q)        
    model<-list(ar = rho1,ma=rho2)
    if (is.null(groups)){
     eps<-filter.arima(model,object$residuals)[-n]
 #print("filter realizado")    
     eps2<-c(eps[(n-maxpq):(n-1)],numeric(nn))
     n<-length(object$residuals)
     x2<-c(object$residuals[(n-maxpq+1):n],numeric(nn))      
     nmax<-maxpq+nn
     for (i in (maxpq+1):nmax)   {
       x2[i] <- sum(model$ar * x2[(i-1):(i-p)]) + ifelse(is.null(model$ma),0,sum(c(1,model$ma)*eps2[i:(i-q)]))    
       }
     x2<-x2[(maxpq+1):nmax]
     ype<-x2
     #yp<-ypx+ype
#print(model) 
# print(p)
#print(q)
#  print("#alternativa a este rollo ")
      ab<-arima(object$residuals, order = c(p,0,q),fixed=coef(object$modelStruct,unconstrained=FALSE)
          ,include.mean=FALSE)
      ab<-predict(ab,n.ahead=n.ahead)

# print("a mano")
# print(ype)
# print("pred arima")
# print(ab)
# print("pred regre")
# print(ypx)
# print("petar")
yp<-as.numeric(ypx)+as.numeric(ab$pred)
# print("petarrr")
# print(ype)
# print(ab)
# print("yp")
# print(yp)
     }
     else{
 #    print("Conditional prediction for each group are not yet done.")
      lengr<-nlevels(gr)
      name.group<-levels(gr)     
      for (i in 1:length(cls1)){
#            cat(" entra grupo ",i)      
            ii.old<-gr.old==name.group[i]
            ng.old<-sum(ii.old)
            ii<-gr==name.group[i]
            ng<-sum(ii)            
            ype2<-numeric(ng)        
            order.ind<-order(ind.old[[i]])
            order.ind.old<-order(ind.old[[i]])            
            res.gr<-object$residuals[ii.old][order.ind.old]
            eps<-filter.arima(model,res.gr)[-ng.old]#el ultimo es NA, por tanto tengo n-1 residuos
   #         print("filter realizado grupo i")    
            eps2<-c(eps[(ng.old-maxpq):(ng.old-1)],numeric(ng))
            #n<-length(object$residuals)
            x2<-c(res.gr[(ng.old-maxpq+1):ng.old],numeric(ng))      
            nmax<-maxpq+ng
            for (j in (maxpq+1):nmax)   {
               x2[j] <- sum(model$ar * x2[(j-1):(j-p)]) + ifelse(is.null(model$ma),0,sum(c(1,model$ma)*eps2[j:(j-q)]))
            }
           x2<-x2[(maxpq+1):nmax]
# print("ype");print(ype)           
# print(ii)
# print(sum(ii))
# print(x2)
# print(order.ind)    
           ype[ii]<-x2#[order.ind]
      }
         yp<-ypx+ype 
         }  
     },
    "corExp"={          
#     else { #gls spatial
     yvalues<-nlme::getCovariate(object$correlation,data=data$df)
#     loci <- gr<-newx$df[,object[["correlation"]][[1]][["index"]]] #buscar en form de gls ~x+y
     rango<-coef(object$modelStruct[[1]],FALSE)  
     sigmasq<- object$sigma^2
# print("rango 1")
# print(rango)     
# print("sigmasq")
# print(sigmasq)   
   cls1<-nlme::getCovariate(object$correlation,data=data$df)  
#   (attributes(res2$modelStruct[[1]])$formula[[2]][[2]])
    gr<-  nlme::getGroups(object$correlation,data=data$df)   
    ind<- nlme::getCovariate(object$correlation,data=data$df)
    ind.old<- nlme::getCovariate(object$correlation,data=object$data$df)
    gr.old<-  nlme::getGroups(object$correlation,data=object$data$df)    
    name.group<-levels(gr) 
#       print(gr);    print("gr")
#      print(ind)    ;    print("ind")
#     print(gr.old);    print("gr.old")
#     print(ind.old)    ;    print("ind.old") 
#   print(name.group)
# print(aemet.sant2$df$longitude)
# print(aemet.sant2$df$latitude)
#  print("lat")
#    print(cls1)     
    if (!is.null(object$correlation))       groups <- nlme::getGroupsFormula(object$correlation)
    else groups <- NULL
#    print(groups)

lon<-paste(nlme::getCovariateFormula(object$correlation)[[2]][[2]]) 
# print("lat00")
lat<-paste(nlme::getCovariateFormula(object$correlation)[[2]][[3]])
#print(coo)
# print(lat)
# print("lat00")
# print(length(ind.old))          
# print(length(gr.old))
#       print("gr.old")
      
      
#hace una nueva prediccion completa de los parametros
# porer si se quiere un rango piferente para cada suboncunto
#     kc1 <- krige.conv(data=object$residuals,coords=coo,locations=newcoo,krige=krige.control(cov.pars=c(.5,1)))
#kc2 <- krige.conv(data=object$residuals,coords=coo,
#locations=newcoo,krige=krige.control(cov.pars=c(sigmasq,rango)))
#     ype<-kc2$predict
#     print(ype)
levgr2<-levels(gr)
levgr<-levels(gr.old)
for (i in 1:length(levgr2)) {
# print(i)

     ii<-gr.old==levgr[i]
     ii2<-gr==levgr[i]     
# print(sum(ii))
# print(sum(ii2))
     coo<-object$data$df[ii,c(lon,lat)]
     newcoo<-    data$df[ii2,c(lon,lat)]
     distxy<-as.matrix(dist(rbind(coo,newcoo),diag=T,upper=T))
n1<-nrow(coo)
n<-n1+nrow(newcoo)                                                       
sig1<-sigmasq*exp(-distxy/rango)
#sig1<-cov.spatial(distxy,cov.model="exponential",cov.pars=c(sigmasq,rango))                        
 #print(sig1[1:24,1:5])
 #print(sig10[1:24,1:5])
 #print("sittttttttttttttttttttttttttttttttt")
sig22<-sig1[1:n1,1:n1] 
sig12<-sig1[(n1+1):n,1:n1]
sig21<-sig1[1:n1,(n1+1):n]
sig11<-sig1[(n1+1):n,(n1+1):n]
 #print(sigmasq)
 # print(rango)
#mu1_2<-   t(sig21)%*%solve(sig22)%*%as.matrix(coo)#valta incluir la parte(X2-mu2)
# print("muuuuuuuuuuuuuuuu000000000")
#  print(dim(sig12))
#  print(dim(sig22))
W0<-sig22
    W <- try(solve(W0),silent=TRUE)
    if (class(W)=="try-error") {
      sv<-svd(W0)
      W<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
      warning("Inverse of sigma computed by SVD")
      }
# print(dim(W))
#ype<- drop(t(sig21)%*%W%*%as.matrix(object$residuals-mean(object$residuals)))
# print("ok1")
bb<-as.matrix(object$residuals[ii])
# print(dim(bb))
aaa<-drop(t(sig21)%*%W%*%as.matrix(object$residuals[ii]))
# print("ok2")
# print(aaa)
# print(sum(ii))
# print(sum(ii2))
ype[ii2]<- drop(t(sig21)%*%W%*%as.matrix(object$residuals[ii]))
}
#ype<- drop(t(sig21)%*%solve(sig22)%*%as.matrix(object$residuals))
#print("muuuuuuuuuuuuuuuu12")
#sig1_2<-diag(sig11-sig12%*%solve(sig22)%*%sig21)
#print("sigmaaaaaaaaaaaaaaa12")
#print(sig1_2)
#print(mu1_2)
#print(cbind(ypx,ype,kc2$predict,ype-kc2$predict,ypx+ype-data$df$yprec))
# print(cbind(ypx,ype,ypx+ype-data$df$yprec))






#print("geoR")
#print(names(kc2))
#    loci<-SpatialPoints(loci)
#    sim <- krige(formula = residual~1, b$df,loci,model=b$fit)
#print(names(sim))    
#    ype<-sim[1]$var1.pred    
#print("gstat")


     })#end switch
#     return(list(yp,ype))

predictor<-drop(yp)
res.var<-object$sr2 
if (missing(pred.var)) pred.var = res.var/weights
if (se.fit || interval != "none") {    
  XX2<-as.matrix(XX[,colnames(object$Vp),drop=FALSE])
   # print("peeeetaaa1")    
  #   print((XX2))
  # print(object$Vp)  
  ip<-rowSums((XX2 %*%object$Vp)*XX2)   
  # print("peeeetaaa2")      
  df<-object$df.residual
  print(ip)
  print(ab)
  se<-sqrt(ip+ab$se^2)#+ype[[2]]                        
  # print(ip)
  # print(ab$se^2)
  yp<-list("pred"=yp,"se"=se)
  
#  if (interval != "none") {
#    tfrac <- qt((1 - level)/2, df)
#hwid <- tfrac * switch(interval, confidence = sqrt(ip),prediction = sqrt(ip + pred.var))    
#predictor <- cbind(predictor, predictor + hwid %o%c(1, -1))
#colnames(predictor) <- c("fit", "lwr", "upr")
#}
}

#else yp<-ypx+ype
     return(yp)
# print("sale2")
}
#     print(names(XX))
#     print(object$formula)
# print(se.fit)
     yp=predict(object,newdata=XX,type=type,se.fit=se.fit,interval=interval,level=level,weights=weights,pred.var=pred.var,df=df,scale=scale,...)
# print(12)
     if (is.null(newx$corStruct)) {
        ype=predict(object$corStruct,se.fit=se.fit,n.ahead=nn)
        if (se.fit) {
          yp[[1]]<-yp[[1]]+ype[[1]]  
          yp[[2]]<-yp[[2]]+ype[[2]]  
          }
        else yp<-yp+ype
     }  
     return(yp)
    #return(predict(object=object,newdata=XX,...))   
#print("ea ea ea")      
#    return(predict(object=object,newdata=XX,type=type,
#    se.fit=se.fit,interval=interval,level=level,weights=weights,pred.var=pred.var,df=df,scale=scale,...))  #cambiarlo    
      }
  else {  
#si pc o pls    
#                                                                                                                         print("ea ea ea") ############  
  for (i in 1:length(vfunc)){
#  if (object$call[[1]]=="fregre.pls")  return(predict(object=object,newdata=XX,type=type,se.fit=se.fit,...)) # poner a yp para despues sumar la parte temp/espacial
  if (object$basis.x[[vfunc[i]]]$type=="pc" |object$basis.x[[vfunc[i]]]$type=="pls") {
#   print("pls pc22")
    a1<-object$coefficients[1]*rep(1,len=nrow(XX))
   object$beta.l[[vfunc[i]]]$data<-matrix(object$beta.l[[vfunc[i]]]$data,nrow=1)
   b1<-inprod.fdata(fdata.cen(newx[[vfunc[i]]],object$mean.list[[vfunc[i]]])[[1]],object$beta.l[[vfunc[i]]])
   yp<-yp+b1
#   yp2<-a1+b1 
#   XX2<-cbind(rep(1,len=nn),XX)
# print("sale pls pc22")
   }
   else{
# print("vuelve fda basis")   
    xcen<-fdata.cen(newx[[vfunc[i]]],object$mean.list[[vfunc[i]]])[[1]]
    x.fd=Data2fd(argvals=tt,y=t(xcen$data),basisobj=object$basis.x[[vfunc[i]]])
    C=t(x.fd$coefs)
    cnames<-colnames(object$JJ[[vfunc[i]]] )
#print(cnames)
    b.est<-matrix(object$coefficients[cnames],ncol=1)
#print( dim(C))
#print(object$JJ[[vfunc[i]]]   )
#print(dim(b.est))
#print(b.est[cnames])
#print(b.est)
    b1<- C%*%object$JJ[[vfunc[i]]]%*%b.est
#print("ha petado")    
#print(yp)
#print(b1)
    yp<-yp+b1  
   }
  } 
  #XX2<-as.matrix(cbind(rep(1,len=nn),XX) )
                        XX2<-as.matrix(XX) #ojo si viene sin intercept hay que descomentar la linea anterior
#   if (se.fit and pc) {
#     se.fit<-sqrt(rowSums((XX2 %*%object$Vp*XX2)))
#     return(list("fit"=yp,"se.fit"=se.fit))
#    } 
  predictor<-drop(yp)
  res.var<-object$sr2 
  if (se.fit || interval != "none") {    
XX2<-as.matrix(XX)
# print("peeeetaaa")    
# print((XX2))
# print(object$Vp)
    ip<-rowSums((XX2 %*%object$Vp)*XX2)   
# print("peeeetaaa222")    
# print(ip)

#    res.var<-sum(fit$residuals^2)/fit$df.residual
       
    df<-object$df.residual
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = sqrt(ip),prediction = sqrt(ip + pred.var))    
        predictor <- cbind(predictor, predictor + hwid %o%c(1, -1))
        colnames(predictor) <- c("fit", "lwr", "upr")
    }
  } 
# print("antes corSTruct")             
if (!is.null(object$corStruct)) {    
# print("entra corSTruct")             
if (names(object$correlation)=="cor.AR"|names(object$correlation)=="cor.ARMA")  {
# print("entra cor.AR cor.ARMA") 
if ((class(object$corStruct[[1]])[1]=="Arima" | class(object$corStruct[[1]])[1]=="ar") & length(object$corStruct)>1)  {
    ype<-NULL
# print("para cada grupo")
#print(object[["correlation"]][[1]][["group"]])
    gr<-newx$df[,object[["correlation"]][[1]][["group"]]]
    lev<-levels(gr) #
    tab<-table(gr)    
    for (j in 1:length(tab)){   
# print(j)    
      lennn<-tab[j]
      previousone<-object$corStruct[[lev[j]]]
#     ind<-gr==lev[j]
      if (lennn!=0) {
#        ype[ind]=predict(object$corStruct[[lev[i]]],object$residuals[ind],se.fit=se.fit,n.ahead=lennn)  
if (class(object$corStruct[[1]])[1]=="Arima")    ype[ind]=predict(object$corStruct[[lev[j]]],se.fit=se.fit,n.ahead=lennn)
#####if (class(object$corStruct[[1]])[1]=="ar")       ype[ind]<-0 # si es diferente viene del arima y no del ar
#print(ype[ind])
#print(object$corStruct[[lev[j]]])
        #ype[ind]=predict(object$corStruc[[lev[i]]],se.fit=se.fit,n.ahead=lennn)           
      }
    } #print("cor Struct AR")  
    }    
# print(names(object$corStruc))              
# print(class(object$corStruc$ar))
if (class(object$corStruct$ar)=="Arima")  ype=predict(object$corStruct$ar,se.fit=se.fit,n.ahead=nn)
if (class(object$corStruct$ar)=="ar")      ype=predict(object$corStruct$ar,object$residuals,se.fit=se.fit,n.ahead=nn)
#ype=NULL# PQ DEBE SER 0predict(object$corStruct$ar,se.fit=se.fit,n.ahead=nn)
#print(ype);print("ype3")
#    ype=predict(object$corStruct$ar,object$residuals,se.fit=se.fit,n.ahead=nn)
#print("ype")    ;print(ype)
if (class(object$corStruct[[1]])[1]=="lm")  {
#    print("cor Struct lm")
    coef.lm<-coef(object$corStruct$lm)
    p<-length(coef.lm)
    lenp<-length(p)
    ype<-NULL
#print("para cada grupo")
#print(object[["correlation"]][[1]][["group"]])
    gr<-newx$df[,object[["correlation"]][[1]][["group"]]]
    if (!is.factor(gr)) gr<-factor(gr)
    lev<-levels(gr) #
    tab<-table(gr)
   for (j in 1:length(tab)){
    
#print("entra en cada grupo para el lm")     
#print(coef.lm)     
      lennn<-tab[j]
      
      previousone<-object$corStruct$lm$res.x[,j]
      #  ahora solo var si estan en el mismo ordenss
      ind<-gr==lev[j]
      e<-rep(NA,len=lennn)
      
#print(ind)  
    
      for (i in 1:lennn) { 
#   print("i")
#   print(i)     
        e[i]<- sum(coef.lm * previousone)
# print(previousone)       
        previousone<-c(e[i],previousone[-p])
#print(e[i])        
#print(previousone)
#        cat("  cor Struct2")
        }
#          print(e)
#        e<-0
#print(ind)
#print(ype)
        ype[ind]<-e         
       }
    }
}       
# 2016_10_31 1.2.4 se comenta esta opcion
#  if (names(object$correlation)=="corVgm")  {
  #    b<- object$corStruct
  # if (b$fit$range[2]<=0) ype<-NULL
  # else{
  # loci <- gr<-newx$df[,object[["correlation"]][[1]][["index"]]]
  # # b<-corVgm(xy=xy,df=df2) 
  # loci<-SpatialPoints(loci)
  ##    gridded(loci) = ~x+y       
  # sim <- krige(formula = residual~1, b$df,loci,model=b$fit)
  # ype<-sim[1]$var1.pred
  # }
#fdalta incluir grupos   
  # }    
#print("corgeoR")
######################################################################################
# se comenta esta opcion par ano tener que cargar la libreria geoR funcion krige.conv  
# buscar fichero previo al 31/10/2016 para utilizar dicho codigo 
# #print(names(object$correlation))  
# if (names(object$correlation)=="cor.Exp")  {
# #print("ook1")
# a<- object$corStruct         
# if (class(a[[1]])[1]=="likGRF")  {
  # #print("ook2")   
  # loci <- gr<-newx$df[,object[["correlation"]][[1]][["index"]]]
# #a<-corgeoR(df2)
# #b<-corVgm(xy=xy,df=df2)     
# #print(loci)        
# #print(names(a))                       
# #kc1 <- krige.conv(a$df, loc=loci,krige=krige.control(cov.pars= a$fit$cov.pars))
# kc2 <- krige.conv(a$df,loc=loci,krige=krige.control(cov.pars=c(a$likfit$sigmasq,a$likfit$phi)))
# #r1<-sum((kc1$predict-sim1$data[51:100])^2)
# #r2<-sum((kc2$predict-sim1$data[51:100])^2)
# #print(kc2)
# #print("residuals espaciales")
# #print(kc1$predict)
# ype<-kc2$predict
# #ype<-sim1$data[(n1+1):n]
# }
# #print("kooo3")          
# } 
  ######################################################################################
  # if (names(object$correlation)=="corCloud")  {
  # #print("ook1")
  # a<- object$corStruct         
  # if (class(a[[1]])[1]=="variomodel")  {
  # #print("ook2")   
  # #a<-corgeoR(df2)
  # #b<-corVgm(xy=xy,df=df2)     
  # #print(loci)        
  # #print(names(a))                       
  # #kc1 <- krige.conv(a$df, loc=loci,krige=krige.control(cov.pars= a$fit$cov.pars))
  # gr2<-newx$df[,object[["correlation"]][[1]][["group"]]]
  # gr<-object$data$df[,object[["correlation"]][[1]][["group"]]]
  # loci <- object$data$df[,object[["correlation"]][[1]][["index"]]]
  # loci2 <- newx$df[,object[["correlation"]][[1]][["index"]]]
  # tab<- table(gr2)
  # lev<-names(tab)
  # #print(gr)
  # #print(gr2)
  # #print(tab)
  # #print(lev)
  # #print("aaaa")
  # ype<-numeric(length(predictor))
  # 
  # for (j in 1:length(tab)) {
  # ind<-gr==lev[j]
  # ind2<-gr2==lev[j]
  # #  kc2 <- krige.conv(a$df[ind,],loc=loci[ind,]krige=krige.control(cov.pars=c(a$likfit$sigmasq,a$likfit$phi)))
  # #     print(loci[ind,])
  # #     print(loci2[ind2,])       
  # #print(a[[1]]$cov.pars)
  # a[[1]]$cov.pars[2]<-a[[1]]$cov.pars[2]*3
  # kc2 <- krige.conv(data=object$residuals[ind],coords=loci[ind,],
  # locations=loci2[ind2,],krige=krige.control(cov.pars=a[[1]]$cov.pars))
  # #
  # #r1<-sum((kc1$predict-sim1$data[51:100])^2)
  # #r2<-sum((kc2$predict-sim1$data[51:100])^2)
  # #print(kc2)
  # #print("residuals espaciales")
  # #print(kc1$predict)
  # ype[ind2]<-kc2$predict
  # #print(kc2$predict)
  # }
  # #ype<-sim1$data[(n1+1):n]
  # }
  # #print("kooo3")          
  # } 
# 31/10/2016fin comentado   names(object$correlation)=="corCloud"  
###################################################################################### 
#else {  print("no correlaciones encontradas")   }
#print(ype)
#print(predictor)
#### fin correlaciones separables      
#print(predictor)
    if (se.fit) {
          predictor<-predictor+ype[[1]]  
          se<-sqrt(ip)+ype[[2]]  
          }
  else predictor<-predictor+ype
  }
  else {   if (se.fit)    se <- sqrt(ip) }
  if (se.fit)   {
        return(list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var)))
        }
  else return(predictor)          
 }    
 }
return(predictor)  
}
#####