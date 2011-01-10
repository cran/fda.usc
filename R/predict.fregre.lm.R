predict.fregre.lm<-function(object,newx=NULL,...){
 if (is.null(object)) stop("No fregre.lm object entered")
 if (is.null(newx)) {
    print("No newx entered")
    newx=object$data # los XX
    }
 data=newx
 basis.x=object$basis.x
 basis.b=object$basis.b
 formula=object$formula
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
if   (attr(tf,"intercept")==0) {
     print("No intecept")
     pf<- paste(pf,-1,sep="")
     }
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
if (length(vnf)>0) {
 first=FALSE
 XX=data.frame(data[[1]][,c(vnf2)])
 names(XX)=vnf2 #data.frame el 1er elemento de la lista
 for ( i in 1:length(vnf)){
# print(paste("no functional variable",vnf[i]))
     if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
if   (attr(tf,"intercept")==0) {
     print("No intecept")
     pf<- paste(pf,-1,sep="")
     }
}
else  first=TRUE
if (length(vfunc)>0)  {
#print(paste("vfunc",vfunc))
        for (i in 1:length(vfunc)) {
#print("class fd")
        	if (class(data[[vfunc[i]]])[1]=="fd"){
        	x.fd=data[[vfunc[i]]]
        	r=x.fd[[2]][[3]]
        	tt = seq(r[1], r[2], len = length(x.fd[[3]]$time))
          basis.x[[vfunc[i]]]=x.fd[[2]]
if (!any(names(basis.b)==vfunc[i])) basis.b[[vfunc[i]]]=x.fd[[2]] #NEW
            Psi = getbasismatrix(tt, basis.x[[vfunc[i]]])
            Theta = getbasismatrix(tt, basis.b[[vfunc[i]]])
            J = t(Psi) %*% Theta
            Z = t(x.fd$coefs) %*% J
            colnames(Z) = paste(vfunc[i], ".",substr(basis.b[[vfunc[i]]]$type,1,3), 1:ncol(Z), sep = "")
            if (first) {
 #             print("no vnf")
              XX=Z
              first=FALSE
              }
            else XX = cbind(XX, Z)
        	for ( i in 1:length(colnames(Z))){if (kterms >= 1)
                  pf <- paste(pf, "+", colnames(Z)[i], sep = "")
                else pf <- paste(pf, colnames(Z)[i], sep = "")
                kterms <- kterms + 1}
        	}
     else {
      if (!is.fdata(data[[vfunc[i]]])) fdataobj=fdata(data[[vfunc[i]]])
      else fdataobj=data[[vfunc[1]]]
      x.fd<-fdataobj[["data"]]
      tt<-fdataobj[["argvals"]]
      rtt<-fdataobj[["rangeval"]]
      if (!any(names(basis.x)==vfunc[i])) {
           np=nrow(x.fd)
           nbasis=ifelse(floor(np/3)>floor((np -4)/2),floor((np-4)/2),floor(np/3))
           basis.x[[vfunc[i]]]=create.bspline.basis(rangeval=rtt,nbasis=nbasis)
       }
				x.fd = Data2fd(argvals = tt, y = t(matrix(x.fd,ncol=length(tt))), basisobj = basis.x[[vfunc[i]]],fdnames=rownames(x.fd))
        		r=x.fd[[2]][[3]]
            Psi = getbasismatrix(tt, basis.x[[vfunc[i]]])
            Theta = getbasismatrix(tt, basis.b[[vfunc[i]]])
            J = t(Psi) %*% Theta
            Z = t(x.fd$coefs) %*% J
            colnames(Z) = paste(vfunc[i], ".",substr(basis.b[[vfunc[i]]]$type,1,3), 1:ncol(Z), sep = "")
            if (first) {
 #              print("no vnf")
              XX=Z
              first=FALSE
              }
            else XX = cbind(XX, Z)
        	  for ( i in 1:length(colnames(Z))){if (kterms >= 1)
                  pf <- paste(pf, "+", colnames(Z)[i], sep = "")
                else pf <- paste(pf, colnames(Z)[i], sep = "")
                kterms <- kterms + 1}
				}
        }
        }
if (is.null(newx)) yp=predict.lm(object,...)
else  {
  if (!is.data.frame(XX)) XX=data.frame(XX)
#  print(names(XX))
  yp=predict.lm(object,XX,...)
}
return(yp)
}


