fregre.lm=function(formula,data,basis.x=NULL,basis.b=NULL,...){
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
name.coef=nam=par.fregre=beta.l=list()
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
 names(XX)=response
}
print(paste("Functional covariate:",vfunc))
if (length(vfunc)>0) {
 mean.list=vs.list=JJ=list()
 bsp1<-bsp2<-TRUE
 for (i in 1:length(vfunc)) {
	if(class(data[[vfunc[i]]])[1]=="fdata"){
      tt<-data[[vfunc[i]]][["argvals"]]
      rtt<-data[[vfunc[i]]][["rangeval"]]
      fdat<-data[[vfunc[i]]];      dat<-data[[vfunc[i]]]$data
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-create.fdata.basis(fdat,l=1:7)
      else   if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)  basis.b[[vfunc[i]]]<-create.fdata.basis(fdat)
      else           if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp2=FALSE
      if (bsp1 & bsp2) {
          if (is.null(rownames(dat)))    rownames(fdat$data)<-1:nrow(dat)
          fdnames=list("time"=tt,"reps"=rownames(fdat[["data"]]),"values"="values")
          xcc<-fdata.cen(data[[vfunc[i]]])
          mean.list[[vfunc[i]]]=xcc[[2]]
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
              }
#          else int<-1:basis.x[[vfunc[i]]]$nbasis
#                  nam[[vfunc[i]]]<-  basis.x[[vfunc[i]]]$names
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
              }
    	    x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]],fdnames=fdnames)

          r=x.fd[[2]][[3]]
#          Theta = getbasismatrix(tt,basis.b[[vfunc[i]]])
#          Psi = getbasismatrix(tt,basis.x[[vfunc[i]]])
#          J = t(Psi) %*% Theta
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
          XX = cbind(XX,Z)
          for ( j in 1:length(colnames(Z))){
           if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
           else pf <- paste(pf, colnames(Z)[j], sep = "")
           kterms <- kterms + 1
           }
        	JJ[[vfunc[i]]]<-J
				}
      else {
        l<-nrow(basis.x[[vfunc[i]]]$basis)
        vs <- t(basis.x[[vfunc[i]]]$basis$data)
        Z<-basis.x[[vfunc[i]]]$x
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",rownames(basis.x[[vfunc[i]]]$basis$data),sep ="")
        XX = cbind(XX,Z)
        vs.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$basis
        mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$mean
        for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }
        #   pf=paste(pf,"-1")
      }
    }
 	else {
 		if(class(data[[vfunc[i]]])[1]=="fd"){
      fdat<-data[[vfunc[i]]]
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-fdat$basis
      else   if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)
         basis.b[[vfunc[i]]]<-create.fdata.basis(fdat,
         l=1:max(5,floor(basis.x[[vfunc[i]]]$nbasis/5)),type.basis=basis.x[[vfunc[i]]]$type,
         rangeval=fdat$basis$rangeval)
      else           if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp2=FALSE
      if (bsp1 & bsp2) {
          r=fdat[[2]][[3]]
#          tt = seq(r[1], r[2], len = length(fdat[[3]]$time))
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
              }
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
              }
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          mean.list[[vfunc[i]]]<-mean.fd(x.fd)
          x.fd<-center.fd(x.fd)
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
          XX = cbind(XX,Z)
          for ( j in 1:length(colnames(Z))){
           if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
           else pf <- paste(pf, colnames(Z)[j], sep = "")
           kterms <- kterms + 1
           }
        	JJ[[vfunc[i]]]<-J
				}
      else {
        l<-ncol(basis.x[[vfunc[i]]]$scores)
        vs <- basis.x[[vfunc[i]]]$harmonics$coefs
        Z<-basis.x[[vfunc[i]]]$scores
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
        XX = cbind(XX,Z)
        vs.list[[vfunc[i]]]=vs
        mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$meanfd
        for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }
         }
    }
   else stop("Please, enter functional covariate")
   }
  }  }
 if (!is.data.frame(XX)) XX=data.frame(XX)
    par.fregre$formula=pf
    par.fregre$data=XX
    z=lm(formula=pf,data=XX,x=TRUE,y=TRUE,...)
#    z$call<-z$call[1:2]
for (i in 1:length(vfunc)) {
 if (bsp1) beta.l[[vfunc[i]]]=fd(z[[1]][name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]])
 else{
	if(class(data[[vfunc[i]]])[1]=="fdata"){
     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
#     beta.est$data<-apply(beta.est$data,2,sum)
     beta.est$data<-colSums(beta.est$data)
     beta.est$names$main<-"beta.est"
     beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
     beta.l[[vfunc[i]]]<-beta.est
     }
 else {
     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*t(vs.list[[vfunc[i]]])
#     beta.est<-apply(beta.est,2,sum)
     beta.est<-colSums(beta.est)
     beta.l[[vfunc[i]]]<-fd(beta.est,basis.x[[vfunc[i]]]$harmonics$basis)
      }
}
}

 z$beta.l=beta.l
 z$formula=pf
 z$mean=mean.list
 z$formula.ini=formula
 z$basis.x=basis.x
 z$basis.b=basis.b
 z$JJ<-JJ
 z$data=z$data
 z$XX=XX
 z$data<-data
 z$vs.list=vs.list   ##### transformarlo en Data2fd(tt,vs.list,basisobj=basisobj)
 class(z)<-c(class(z),"fregre.lm")
 z
}


