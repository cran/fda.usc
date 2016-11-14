
###############################################################################
###############################################################################
fregre.igls=function(formula,data,basis.x=NULL,basis.b=NULL,correlation,maxit=100,
rn,lambda,weights=rep(1,n),control,...){
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
# vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 name.coef=nam=par.fregre=beta.l=list()
 kterms=1
 n<-length(data[["df"]][,response])
 XX=data.frame(data[["df"]][,c(response)],weights)
 namxx=names(XX)=c(response,"weights")
 aa<-NULL
 if (length(vnf)>0) {          
#print(paste("Non functional covariate:",vnf))
# XX=data.frame(XX,data[["df"]][,c(vnf)])    # si es facotr poner k-1 variables dummies
# names(XX)=c(namxx,vnf)
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
#  contrasts = NULL
    mf<-as.data.frame(model.matrix(formula(pf),data$df))
    vnf2<-names(mf)[-1]
    pf <- rf <- paste(response, "~", sep = "")
    for ( i in 1:length(vnf2))  pf<-paste(pf, "+", vnf2[i], sep = "")
    cname<-names(XX)
    XX <- data.frame(XX,mf)
    names(XX)<-c(cname,names(mf))
}
else {
    XX <- data.frame(XX,model.matrix(formula(paste(pf, "1")),data$df))
    names(XX)[3]<-"(Intercept)"
}
#print(paste("Functional covariate:",vfunc))
if (missing(rn))    {    rn0=FALSE;                    rn=list()}
else rn0<-TRUE
if (missing(lambda))    {    lambda0=FALSE;                    lambda=list()}
else lambda0<-TRUE
mat<-rep(0,len=ncol(XX)-2)
imat2<-ncol(XX)-2
mat2<-diag(0,nrow=imat2)
if (length(vfunc)>0) {
 mean.list=vs.list=JJ=list()
 bsp1<-bsp2<-TRUE
 for (i in 1:length(vfunc)) {
	if(class(data[[vfunc[i]]])[1]=="fdata"){
      tt<-data[[vfunc[i]]][["argvals"]]
      rtt<-data[[vfunc[i]]][["rangeval"]]
      np<-length(tt)
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
         nbasisb<- basis.b[[vfunc[i]]]$nbasis
         R=diag(0,ncol= nbasisb,nrow=nbasisb)
         R<-eval.penalty(basis.b[[vfunc[i]]],Lfdobj = int2Lfd(2))
#        mat2<-diag(0,nrow=imat2+nbasisb)
        lenm<-ncol(mat2)
        if (rn0) stop("Ridge regressions is only implemented for functional principal component basis")
        MM<-matrix(0,nrow=lenm,ncol=nbasisb)       
        MM2<-matrix(0,nrow=nbasisb,ncol=imat2+nbasisb)
        mat2<-cbind(mat2,MM)         
        mat2<-rbind(mat2,MM2)               
        if (!is.null(lambda[[vfunc[i]]]))  {
          mat2[(imat2+1):(imat2+nbasisb),(imat2+1):(imat2+nbasisb)]<-lambda[[vfunc[i]]]*R
          
        } 
        imat2<-imat2+nbasisb    
			}
      else {
        basis<-basis.x[[vfunc[i]]]
        l<-basis$l
        vs <- t(basis$basis$data)           
        basis$x<-basis$x[,l,drop=FALSE]
        Z<-basis$x
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",rownames(basis$basis$data),sep ="")
        XX = cbind(XX,Z)
        vs.list[[vfunc[i]]]=basis$basis
        mean.list[[vfunc[i]]]=basis$mean
        for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }
        if (!is.null(rn[[vfunc[i]]]))  {
           mat<-c(mat,rep(rn[[vfunc[i]]],len=length(l)))
           }     
        else          mat<-c(mat,rep(0,len=length(l)))
#        mat2<-diag(0,nrow=imat2+length(l))
        lenl<-length(l)
        lenm<-ncol(mat2)
        MM<-matrix(0,nrow=lenm,ncol=lenl)       
        MM2<-matrix(0,nrow=lenl,ncol=imat2+lenl)                
        mat2<-cbind(mat2,MM)         
        mat2<-rbind(mat2,MM2)         
        if (!is.null(lambda[[vfunc[i]]]))  {
#        mat2<-diag(0,nrow=imat2+length(l))
        R<-P.penalty(1:length(l))
        mat2[(imat2+1):(imat2+lenl),(imat2+1):(imat2+lenl)]<-lambda[[vfunc[i]]]*R
        imat2<-imat2+lenl
        }                    
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
#print("acaba creacion de XX")  
 if (!is.data.frame(XX)) XX=data.frame(XX)
    par.fregre$formula=pf
    par.fregre$data=XX
    y<-XX[,1]    
    scores<-as.matrix(XX[,-(1:2)])     
#    scores<-as.matrix(XX)     
    W<-diag(weights) 
    error<-FALSE
    if (missing(correlation)) correlation0<-FALSE
    else {
      if (is.matrix(correlation)) error=FALSE
      else                        error=TRUE
      correlation0<-TRUE
      }
    if (!rn0&!lambda0&!correlation0) {
      W0<-FALSE
      z=lm(formula=pf,data=XX,...) 
      e<-z$residuals
      S<-solve(t(scores)%*%W%*%scores)          
#       S<-diag(coef(summary(z))[,2])
      class(z)<-c(class(z),"fregre.lm")
       mat0<-diag(0,ncol(scores))
    }      
    else {
#       S<-solve(t(scores)%*%W%*%scores)
       mat0<-diag(0,ncol(scores))
       if (lambda0) {
#             print(" mat  2")
#             print(mat2)    
             mat0<-mat2       #incluir pesos solve(W) 
       }             
       if (rn0) {
#             print(" mat  ")
#             print(mat)    
             mat0<-diag(mat)
       }      
       if (rn0 & lambda0) warning("Only ridge penalization is done by rn argument (lambda argument is ignored)")               
#      ddd<-t(scores)%*%W  
#      S<-solve(t(scores)%*%W%*%scores+mat+mat2)       #incluir pesos solve(W)
      W<-diag(weights)%*%W
      S<-solve(t(scores)%*%W%*%scores+mat0)
      Cinv<-S%*%t(scores)%*%W                    #incluir pesos W repetri proceso hasta que no cambie la prediccion   
      ycen = y - mean(y)
      coefs<-Cinv%*%XX[,1]
      z<-list()      
      z$fitted.values<-yp<-drop(scores%*%coefs)
      e<-z$residuals<-XX[,1]- z$fitted.values      
################################################################################
################################################################################
 it<-1
 eps=.001
 err2=sqrt(sum(coefs^2))
 MM<-matrix( 1:n,ncol=1)
 corStruct<-list()
while (error) {   
 W<-W0<-diag(n)
#print(dim(W)) 
 #    aa<-ar(e,order.max=8)
#    if (aa$order>0) {
#     W0<-toeplitz(ARMAacf(aa$ar,lag.max=n-1))    }
if (is.list(correlation)) {
  ncor<-length(correlation)    
  name.cor<-names(correlation) 
  }
else {
ncor<-1
name.cor<-correlation[[1]]
correlation=list(correlation=list())
names(correlation)<-name.cor
}
#name.cor<-names(correlation) 
wei0<-wei<-diag(weights) 
ee<-e
for (icor in 1:ncor) { 
if (!is.list(correlation[[icor]])) {
      name.cor<-correlation[[icor]]
      correlation[[icor]]=list()
      names(correlation)[icor]<-name.cor
      } 
# print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa111")
# print(names(correlation[[icor]]))            
if (is.null(correlation[[icor]][["group"]])) {
# print("no hay grupo 0")
     n<-length(e)
     par.cor<-switch(name.cor[icor],
             "cor.ARMA"={
          cor.name<-"cor.ARMA"
         if (is.null(correlation[[icor]][["index"]])) index<-1:n
         else {
           index.df<-correlation[[icor]][["index"]]
           index<-data[["df"]][,index.df]
           if (is.vector(index) )       e<-e[order(index)]
           else {
             par.ord<-list()
             for (i in 1:ncol(index)) par.ord[[i]]<-index[,i]
             index<-as.numeric(do.call("order",par.ord))
             e<-e[index]
            }
            }
         list("x"=e[index])
         },
        "cor.AR"={
          cor.name<-"cor.AR"
         if (is.null(correlation[[icor]][["index"]])) index<-1:n
         else {
           index.df<-correlation[[icor]][["index"]]
           index<-data[["df"]][,index.df]
           if (is.vector(index) )       e<-e[order(index)]
           else {
             par.ord<-list()
             for (i in 1:ncol(index)) par.ord[[i]]<-index[,i]
             index<-as.numeric(do.call("order",par.ord))
             e<-e[index]
            }
            }
         list("x"=e[index])
         },
         "corExpo"={     
           index.df<-correlation[[icor]][["index"]]
            dxy<-data[["df"]][,index.df]
#           dxy<-correlation[[icor]]$index
           index<-1:length(e)       
           cor.name<-"corExpo"
           list("xy"=dxy)}, 
         "corVgm"={     
           index.df<-correlation[[icor]][["index"]]
           xy<-data[["df"]][,index.df]
#           if (i==1) 
           dxy0<- as.matrix(dist(xy,diag =TRUE, upper = TRUE))
#           dxy<-correlation[[icor]]$index
           index<-1:n       
           cor.name<-"corVgm"                        
           residual<-e[index]
           list("dxy"=dxy0,"xy"=xy,"df"=data.frame(residual))},   
         "cor.Exp"={     
           index.df<-correlation[[icor]][["index"]]
           xy<-data[["df"]][,index.df]
#           if (i==1) 
        #   dxy<-dxy0<- as.matrix(dist(xy,diag =TRUE, upper = TRUE))
#           dxy<-correlation[[icor]]$index
           index<-1:n       
           cor.name<-"cor.Exp"                        
           residual<-e[index]
# print(index.df)           
# print("index.df")
           dff<-data.frame(data[["df"]][,index.df],residual)
           names(dff)<-c(index.df,"residual")
#        print(names(dff))
# print("corgeoR")       
#           list("df"=dff,"dxy"=dxy)},                                 
           list("df"=dff)},                                 
         "corUnstruc"={
            cor.name<-"corUnstruc"
            list(correlation[[icor]][["index"]])
            }
           ) 
      e<-e[index]  
#     print("NO grupo0")
#     print((par.cor))
#     print(cor.name)
#     print(names(par.cor))
#     print("aaaaaaaaaaaaaa")           
     if (!missing(control)) {
       npar<-length(par.cor)
       for (i in 1:length(control)) {
          par.cor[npar+i]<-control[i]
          names(par.cor)[npar+i]<-names(control)[i]
     }
     }
#      print(names(par.cor))
#      print("NO grupo0 names par.cor ****************")     
#      print(par.cor)
     aa<-do.call(cor.name,par.cor)
     W0<-aa[[1]] 
     corStruct<-aa[-1]
#          print(names(aa))
 #    print(names(corStruct))
 #   print(W0[1:3,1:3])
#    print("W0[1:3,1:3]")
     
#print(W0)
# print("fin no grupo")     
#     W0<-W0[index,index]
}
else {
# print(" si hay grupooooo")    
n<-nrow(data$df)
W00<-matrix(0,ncol=n,nrow=n)
#print("si hay grupo/s")
name.group<-correlation[[icor]]$group
#print(name.group)
ncomb<-length(name.group)
if (ncomb>1) {
gr<-(data[["df"]][,name.group[1]])
for (ngr in 2:ncomb)
  gr<-paste(gr,(data[["df"]][,name.group[ngr]]),sep="")
#print(gr)
#print("factor")
gr<-factor(gr,levels=unique(gr))
#print(gr)
#print("factor2")
}
else gr<-data[["df"]][[name.group]]
# print("gr gr gr gr GR GR Gr")
if (!is.factor(gr)) gr<-factor(gr)
  lev<-levels(gr)
  nlev<-length(lev)
#  ee<-e
    par.cor<-switch(name.cor[icor],   
        "cor.ARMA"={
#         print("se llama  a corARMA con grupos")
#          print(correlation[[icor]])
#         if (is.null(correlation[[icor]][["index"]])) index<-1:n
#         else {
#           print(correlation[[icor]]$index)
           index.df<-correlation[[icor]][["index"]]
           index<-data[["df"]][,index.df]
#           print(index.df)
#           print(index)
#           print(is.vector(index))
#                       print(class(index))
#                       print("ok")
e2<-NULL
for (i in 1:nlev) {
 jj<-gr==lev[i]
 e0<-e[jj]
 e2<-cbind(e2,e0[order(index[jj])])
 }       

#           if (is.vector(index) )       e<-e[order(index)]
#           else {
#             par.ord<-list()
#             for (i in 1:ncol(index)) par.ord[[i]]<-index[,i]
#             index<-as.numeric(do.call("order",par.ord))
#             e<-e[order(index),]
#            }
#          }
         cor.name<-"cor.ARMA"
 #        print(e)
         list("x"=e2)                 
# arrreglar para pasar una matriz !!!              
         }, 
              "corCloud"={     
           index.df<-correlation[[icor]][["index"]]
           index<-data[["df"]][,index.df]           
           gr<-data[["df"]][,correlation[[icor]][["group"]]]
           xy<-data[["df"]][,index.df]
#           dxy<-dxy0<- as.matrix(dist(xy,diag =TRUE, upper = TRUE)) # solo calcular 1 vez
#           index<-1:n       
           cor.name<-"corCloud"                        
#           residual<-e[index]
            residual<-e
#print(index.df)           
# print("index.df")
           dff<-data.frame(data[["df"]][,index.df],residual)
#        print(names(dff))           
           names(dff)<-c(index.df,"residual")
#        print(names(dff))
#        print(class(gr))
#  print("corgeoR")       
           list("df"=dff,"gr"=gr)
           },                  
     "cor.Exp"={     
           index.df<-correlation[[icor]][["index"]]
           index<-data[["df"]][,index.df]           
           gr<-data[["df"]][,correlation[[icor]][["group"]]]
           xy<-data[["df"]][,index.df]
#           dxy<-dxy0<- as.matrix(dist(xy,diag =TRUE, upper = TRUE)) # solo calcular 1 vez
#           index<-1:n       
           cor.name<-"cor.Exp"                        
#           residual<-e[index]
            residual<-e
#print(index.df)           
# print("index.df")
           dff<-data.frame(data[["df"]][,index.df],residual)
#        print(names(dff))           
           names(dff)<-c(index.df,"residual")
#        print(names(dff))
#        print(class(gr))
#  print("corgeoR")       
           list("df"=dff,"gr"=gr)
           },                                        
         "corUnstruc"={
            cor.name<-"corUnstruc"  
            index<-1:nrow(dff)
#            print(correlation[[icor]][["index"]])
            par.cor<-list(correlation[[icor]][["index"]])
            }) 
# print("fin Switch")            
#print(length(ee))
#  print("AAAA2")
#      print(cor.name)
# print(!missing(control)) 
     if (!missing(control)) {
       npar<-length(par.cor)
# print(names(par.cor))       
# print(names(control))    
       for (k in 1:length(control)) {
#       print(control[k])
#      print(k)
#       print(names(control)[k])
          par.cor[npar+k]<-control[k]
          names(par.cor)[npar+k]<-names(control)[k]
     }
     }     
#      print(names(par.cor))
     aa<-do.call(cor.name,par.cor)
     W0<-aa[[1]]  
    corStruct<-aa[-1]
# print(dim(W0))     
#print(length(e))
# print("AAAA3")
# print(dim(W0))
# print("bbbbbbbbbbbbbbbbbb")  
# for (i in 1:nlev) {
# jj<-gr==lev[i]
# e0<-e[jj]
# e2<-cbind(e2,e0[order(index[jj])])
# }     
 for (j in 1:nlev) {
     n<-table(gr)[j]
     ilev<-gr==lev[j]
     e<-ee[ilev] 
#print(cbind(e,index))            
#     e<-e[index]
#     print("si grupo0")
#      print(names(par.cor))
#  print("ha petado  0")    
    if (name.cor[icor]=="cor.Exp" | name.cor[icor]=="corCloud" )    W00[ilev,ilev]<-W0
    else            W00[ilev,ilev]<-W0[order(index[ilev]),order(index[ilev])]     
#    else     W00[ilev,ilev]]<-W0     

  }
#  print("AAA444")
W0<-W00

}
# print(dim(wei0))
# print(dim(W0))    
wei0<-wei0%*%W0
} #fin for ncor
W0<-wei0  
#     W0<-wei%*%W0
#print("antes de petar")
#print(dim(W0))
#     print(W0[1:3,1:3])
#         print("W0[1:3,1:3]")
#     W<-solve(W0)      
     W <- try(solve(W0),silent=TRUE)        
     if (class(W)=="try-error") {
     sv<-svd(W0)
     W<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
#     warning("Inverse of sigma computed by SVD")
     }
#     print(W[1:3,1:3])
# print("despues de peta") 
     
    mat0[1,1]<-0
#print(mat0[1:3,1:3])     
#print(W0[1:3,1:3])     
#print(W[1:3,1:3])
    W2<-t(scores)%*%W%*%scores+mat0
#    S<-solve(W2)              
     S <- try(solve(W2),silent=TRUE)
     if (class(S)=="try-error") {
     sv<-svd(W2)        
     S<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
#     warning("Inverse of sigma computed by SVD")
     }    
#print("aaabbb")    
    Cinv<-S%*%t(scores)%*%W                    #incluir pesos W repetir proceso hasta que no cambie la prediccion   
    coefs<-Cinv%*%matrix(y,ncol=1)
    yp<-drop(scores%*%coefs)    
    coefs<-drop(coefs)
    
#    cnames<-names(XX[,-(1:2)])    
#    names(coefs)<-c("Intercept",cnames)
    e<-y-yp         
#    err3<-sum(e^2) 
    err3<-coefs
#    err4=sqrt(sum((coefs-err2)^2))
    err4<-max(abs((err3-err2)/err2))  #e<.01
#print(c(it,err2,err3,err4))
       if (err4<eps) error<-FALSE
       else {
        if (it==maxit) error<-FALSE
        else it<-it+1   
        err2<-err3
        }                                                                                                #cat(" it it  it it")        ;print(it)        
 }  
################################################################################           
#print("aaabbb2")   
      H<-scores%*%Cinv
      df<-traza(H)
      coefs<-drop(coefs)
#      cnames<-names(XX[,-(1:2)])
#      names(coefs)<-c("Intercept",cnames)  
      z$coefficients<-coefs
      z$mean.list<-mean.list
      z$df.residual<-n-df
      z$H<-H
      z$r2 <- 1 - sum(z$residuals^2)/sum(ycen^2)       
#print("aaabbb3")         
      if  (class(basis.x[[vfunc[1]]])=="basisfd") {
        z$call[[1]] = "fregre.basis"
#        z$lambda<-rn
        }
       else {
        if  (basis.x[[vfunc[1]]]$type=="pc")  z$call[[1]] = "fregre.pc"
        if  (basis.x[[vfunc[1]]]$type=="pls")  z$call[[1]] = "fregre.pls"        
        }             
    class(z)<-c("fregre.fd","fregre.lm")
    rdf<-n-df
    sr2 <- sum(e^2)/ rdf
    r2 <- 1 - sum(e^2)/sum(ycen^2)
#print("aaabbb34")      
    r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
    GCV <- sum(e^2)/(n - df)^2
#print("aaabbb4")     
        z$terms<-terms
        z$residuals <- drop(e)
        z$fitted.values <- yp
        z$y <- y
        
        z$rank <- df
        z$df.residual <-  rdf
        Z=cbind(rep(1,len=n),Z)
        colnames(Z)[1] = "(Intercept)"
        std.error = sqrt(diag(S) *sr2)
        t.value = coefs/std.error
#print("aaabbb5")         
        p.value = 2 * pt(abs(t.value), n - df, lower.tail = FALSE)
#print("aaabbb6")         
        coefficients <- cbind(coefs, std.error, t.value, p.value)
#print("aaabbb7")         
        colnames(coefficients) <- c("Estimate", "Std. Error",
            "t value", "Pr(>|t|)")
        z$coefs<-coefficients    
        class(z) <- "lm"
        z$it<-it
        # z$ar<-aa
# z$W<-W
# z$W0<-W0
# z$it<-it
}       
#    z$call<-z$call[1:2]
for (i in 1:length(vfunc)) {
#print(bsp1)
 if (bsp1) beta.l[[vfunc[i]]]=fd(z$coefficients[name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]])
 else{
	if(class(data[[vfunc[i]]])[1]=="fdata"){
#     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
     beta.est$data<-colSums(beta.est$data) 
     beta.est$names$main<-"beta.est"
     beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
     beta.est$names$main<-"beta.est"
     beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
           if  (basis.x[[vfunc[i]]]$type=="pls") {
             if (basis.x[[vfunc[i]]]$norm)  {
              sd.X <- sqrt(apply(data[[vfunc[i]]]$data, 2, var))
              beta.est$data<-  beta.est$data/sd.X
             }      
            }   
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
# print("saleee4")  
 z$sr2<-sum(e^2)/z$df.residual
 z$Vp=z$sr2*S
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
 z$fdataobj<-data[[vfunc[1]]]
 if (correlation0) {
  rn0<-TRUE
  z$corStruct<-corStruct
  z$it<-it
  }
 z$rn<-rn0
 z$lambda<-lambda0
 z$W<-W
 z$W0<-W0       
# z$M<-mat0 #penalty matrix (ridge or curvature penalty)
 z$correlation<-correlation 
 z$correl<-correlation0 
 z$vs.list=vs.list   
 class(z)<-c("lm","fregre.lm")
 z
}
