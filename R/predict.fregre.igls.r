# ojo predict.fregre.fd o predict.fregre.lm
# que pinta el df=df en el call

################################################################################
predict.fregre.igls<-function(object,new.fdataobj=NULL,data,se.fit=FALSE,scale = NULL,df=df,
    interval = "none", level = 0.95,weights = 1, pred.var = res.var/weights,...){
 if (is.null(object)) stop("No object entered")
 if (is.null(new.fdataobj)) { yp<-object$fitted.values
#    yp=predict(object,type=type,se.fit=se.fit,interval=interval,level=level,weights=weights,pred.var=pred.var,df=df,scale=scale,...)    
    print("No newx entered")
    return(yp)
    }
 else {
# print(object$call[[1]])
#object2<-object
#if (object$call[[1]]=="fregre.npgls") object2$call[[1]]<-"fregre.npgls"
#print(object$call[[1]])                 
#predict.fregre.fd<-function(object,new.fdataobj=NULL,se.fit=FALSE,scale = NULL,df=df,
#    interval = "none", level = 0.95,weights = 1, pred.var = res.var/weights, ...){ 
yp<-predict.fregre.lm(object,new.fdataobj,se.fit=se.fit,scale=scale,df=df,interval=interval,level=level,weights=weights,
                      
                      
pred.var=pred.var,...)                 
predictor<-yp
#print(yp)
# print(object$basis.x[[vfunc[i]]])
# print(object$basis.x[[vfunc[i]]]$type)
#print(vfunc)     
#print("2")       
       
                                    #print(object$rn)
# if (!missing(data)) 
#else    nn<-nrow(new.fdataobj) 
 nn <-length(yp)
 #print(nn)
 
 n<-length(object$residuals) 
  predictor<-drop(yp)
  res.var<-object$sr2 
  
if (!is.null(object$corStruct)) {    
#  print("entra corSTruct")             
if (names(object$correlation)=="cor.AR"|names(object$correlation)=="cor.ARMA")  {
#  print("entra cor.AR cor.ARMA") 
if ((class(object$corStruct[[1]])[1]=="Arima" | class(object$corStruct[[1]])[1]=="ar") & length(object$corStruct)>1)  {
    ype<-NULL
#  print("para cada grupo")
#print(object[["correlation"]][[1]][["group"]])
    gr<-data[,object[["correlation"]][[1]][["group"]]]
    lev<-levels(gr) #
    tab<-table(gr)    
   
    for (j in 1:length(tab)){   
#  print(j)    
      lennn<-tab[j]
      previousone<-object$corStruct[[lev[j]]]
      ind<-gr==lev[j]
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
#  print("cor Struct AR")  
if (class(object$corStruct$ar)=="Arima")  ype=predict(object$corStruct$ar,se.fit=se.fit,n.ahead=nn)
#print("petooooo")
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
# print("para cada grupo")
#print(object[["correlation"]][[1]][["group"]])
    gr<-data[,object[["correlation"]][[1]][["group"]]]
    if (!is.factor(gr)) gr<-factor(gr)
    lev<-levels(gr) #
    tab<-table(gr)
   for (j in 1:length(tab)){
    
#print("entra en cada grupo para el lm")     
#print(coef.lm)     
      lennn<-tab[j]
      
      previousone<-object$corStruct$lm$res.x[,j]  #falta coger el del grupo correspondiente
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
#          print(yp)
#        e<-0
#print(ind)
#print(ype)
        ype[ind]<-e         
       }
    }
} 
  predictor<-predictor+ype
#    if (se.fit) {
#          predictor<-predictor+ype[[1]]  
#          se<-sqrt(ip)+ype[[2]]  # no esta definico 'ip' previamente
#          }
#  else predictor<-predictor+ype
  }
 # else {   if (se.fit)    se <- sqrt(ip) }
  #if (se.fit)   {
#        return(list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var)))
#        }
  #else return(predictor)
   return(predictor)          
 }    

return(predictor)  
}
