################################################################################
predict.fregre.kgam=function(object,newx=NULL,type="response",...){
#if type="response" solo p[,"mu"]
#	nvars=ncol(object$effects)-1 #y si no hay Intercept
	namesx<-names(object$result)
	nvars=length(namesx)
  nr=nrow(newx[[namesx[1]]])
	pr=matrix(NA,nrow=nr,ncol=nvars+3)
	colnames(pr)=c(colnames(object$effects),"eta","mu")
	pr[,"Intercept"]=rep(object$effects[1,"Intercept"],nr)
#	for (i in 1:nvars){pr[,i]=predict(object$result[[i]],xlistnew[[i]])-mean(object$result[[i]]$fitted.values)}
	for (i in 1:nvars){
   pr[,i]=predict(object$result[[namesx[i]]],newx[[namesx[i]]])
   }

   if (nr==1) pr[,"eta"]<-sum(pr[,1:(nvars+1)])
	 else pr[,"eta"]=apply(pr[,1:(nvars+1)],1,sum)
#   pr[,"mu"]=object$family$linkinv(pr[,"eta"])
  pr<-switch(type,"response"=object$family$linkinv(pr[,"eta"]),"link"=pr[,"eta"])
	return(pr)
}
