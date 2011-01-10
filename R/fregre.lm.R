fregre.lm=function(formula,data,basis.x=NULL,basis.b=NULL,w=NULL,...){
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
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
if (length(vnf)>0) {
 XX=data[[1]][,c(response,vnf2)] #data.frame el 1er elemento de la lista
 for ( i in 1:length(vnf)){
print(paste("Non functional covariate:",vnf[i]))
     if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
#     else pf <- paste(pf, terms[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
if   (attr(tf,"intercept")==0) {
#     print("No intecept")
     pf<- paste(pf,-1,sep="")
     }
}
else {
 XX=data.frame(data[[1]][,response])
names(XX)=response  }
print(paste("Functional covariate:",vfunc))
if (length(vfunc)>0) {
 for (i in 1:length(vfunc)) {
 	if (class(data[[vfunc[i]]])[1]=="fd"){
#     print(paste("fd object ",vfunc[i]))
   	 x.fd=data[[vfunc[i]]]
   	 r=x.fd[[2]][[3]]
   	 tt = seq(r[1], r[2], len = length(x.fd[[3]]$time))
     basis.x[[vfunc[i]]]=x.fd[[2]]
     if (!any(names(basis.b)==vfunc[i])) basis.b[[vfunc[i]]]=x.fd[[2]] #NEW
     Psi = getbasismatrix(tt, basis.x[[vfunc[i]]])
     Theta = getbasismatrix(tt, basis.b[[vfunc[i]]])
     J = t(Psi) %*% Theta
     x2=t(x.fd$coefs)
     xmean=apply(x2,2,mean)
     xcen=sweep(x2,2,xmean,FUN="-")
     class(xcen)="matrix"
 #    Z = xcen %*% J
     Z = x2 %*% J #parte nueva
     colnames(Z) = paste(vfunc[i], ".fd", 1:ncol(Z), sep = "")
     XX = cbind(XX, Z)
     for ( i in 1:length(colnames(Z))){
       if (kterms >= 1)   pf <- paste(pf, "+", colnames(Z)[i], sep = "")
       else pf <- paste(pf, colnames(Z)[i], sep = "")
       kterms <- kterms + 1}
     }
	if(class(data[[vfunc[i]]])[1]=="fdata"){
      #tt = attr(data[[vfunc[i]]],"argvals")
      #x<-data[[vfunc[i]]][["data"]]
      tt<-data[[vfunc[i]]][["argvals"]]
      rtt<-data[[vfunc[i]]][["rangeval"]]
      if (!any(names(basis.x)==vfunc[i])) {
        np=nrow(data[[vfunc[i]]])
        nbasis=ifelse(floor(np/3)>floor((np -4)/2),floor((np-4)/2),floor(np/3))
        basis.x[[vfunc[i]]]=create.bspline.basis(rangeval=range(tt),nbasis=nbasis)
      }
fdnames=list("time"=tt,"reps"=rownames(data[[vfunc[i]]][["data"]]),"values"="values")
#       fdnames$time=tt
#       fdnames$reps=rownames(data[[vfunc[i]]])
#      	xmean=apply(x,2,mean)
#       xc=matrix(data[[vfunc[i]]],ncol=length(tt)) ####
#       xcen=sweep(xc,2,xmean,FUN="-")
#       class(xcen)="matrix"
#				x.fd = Data2fd(argvals = tt, y = t(xcen), basisobj = basis.x[[vfunc[i]]],fdnames=fdnames)
	x.fd = Data2fd(argvals = tt, y = t(matrix(data[[vfunc[i]]][["data"]],ncol=length(tt))), basisobj = basis.x[[vfunc[i]]],fdnames=fdnames)
       	r=x.fd[[2]][[3]]
        Psi = getbasismatrix(tt,basis.x[[vfunc[i]]])
        Theta = getbasismatrix(tt,basis.b[[vfunc[i]]])
        J = t(Psi) %*% Theta
        x2=t(x.fd$coefs)
      	xmean=apply(x2,2,mean)
        xcen=sweep(x2,2,xmean,FUN="-")
        class(xcen)="matrix"
#        Z = xcen %*% J
        Z = x2 %*% J #parte nueva
        colnames(Z) = paste(vfunc[i], ".",substr(basis.b[[vfunc[i]]]$type,1,3), 1:ncol(Z), sep = "")
        XX = cbind(XX, Z)
        for ( i in 1:length(colnames(Z))){
          if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[i], sep = "")
          else pf <- paste(pf, colnames(Z)[i], sep = "")
          kterms <- kterms + 1}
				}      }       }
  if (!is.data.frame(XX)) XX=data.frame(XX)
  if (is.null(w))	z=lm(formula=pf,data=XX,x=TRUE,y=TRUE,...)
  else  	z=lm(formula=pf,data=XX,x=TRUE,y=TRUE,weights=w) #no run
if (length(vfunc)>0) {
 for (i in 1:length(vfunc)) {
 if ((class(data[[vfunc[i]]])[1]=="fd")||(class(data[[vfunc[i]]])[1]=="fdata")){
    name.coef=paste(vfunc[i],".",substr(basis.b[[vfunc[i]]]$type,1,3),1:basis.b[[vfunc[i]]][[4]],sep="")
  	beta.l[[vfunc[i]]]=fd(z[[1]][name.coef],basis.b[[vfunc[i]]])
            	}}}
 z$beta.l=beta.l
 z$formula=formula
 z$basis.x=basis.x
 z$basis.b=basis.b
 z$data=z$data
 z$XX=XX
 z$pf=pf
 z$Z=Z
 z$x.fd=x.fd
 z
}





