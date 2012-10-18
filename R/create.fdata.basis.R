#######################
#######################
create.fdata.basis<-function(fdataobj,l=1:5,maxl=max(l),type.basis="bspline",
rangeval=fdataobj$rangeval,class.out="fd"){
      aa1 <- paste("create.",type.basis,".basis", sep = "")
      as <- list()
      as$rangeval <- rangeval
      as$nbasis <-maxl
      as$dropind<-setdiff(1:maxl,l)
      basis=do.call(aa1,as)
      if (class.out=="fdata") {
          nam<-basis$names[intersect(1:maxl,l)]
          basis=fdata(t(eval.basis(fdataobj$argvals,basis)),fdataobj$argvals,fdataobj$rangeval)
          rownames(basis$data)<-nam
          basis$type<-type.basis
          basis$nbasis<-maxl
          basis$dropind<-as$dropind
         }
         basis
}
#######################
#######################
#create.pc.basis<-function(fdataobj,l=1:5,norm=TRUE){
#     pc<-fdata2pc(fdataobj,norm=norm,ncomp=max(l))
#     basis=pc$rotation[l,]
#     rownames(basis$data)<-paste("PC",l,sep="")
#no#     if (length(l)==1)   x=as.matrix(pc$x[,l],ncol=1)
#no#     else    x=pc$x[,l]
#      out<-list("basis"=basis,"x"=pc$x,"mean"=pc$mean,
#      "fdataobj.cen"=pc$"fdataobj.cen","fdataobj"=pc$fdataobj,
#      "lambdas"=pc$lambdas,"l"=l,"type"="pc")
#      class(out)<-"fdata.comp"
#return(out)
#}
#######################
#######################
create.pc.basis<-function(fdataobj,l=1:5,norm=TRUE,basis=NULL,rn=0,...){
 tt<-fdataobj$argvals
 rtt<-fdataobj$rangeval
 dropind=NULL
 pc<-fdata2pc(fdataobj,norm=norm,ncomp=max(l),...)
 vs<-pc$rotation$data
# if (is.vector(pc$u)) pc$u<-matrix(pc$u,ncol=1)
 if (length(l)==1) {
   if (is.matrix(pc$u)) {
     pc$u<-matrix(pc$u[,l],ncol=1) 
     vs<-matrix(vs[l,],nrow=1) 
    pc.fdata<-pc$u%*%matrix(pc$d[l])%*%vs     
     }
   else { 
        print(dim(pc$u))
        pc$u<-matrix(pc$u,ncol=1) 
        pc$rotation$data<-matrix(vs[l,],ncol=1)
        pc.fdata<-pc$u%*%(pc$d)%*%vs        
        }

  }
 else  pc.fdata<-pc$u[,l]%*%diag(pc$d[l])%*%vs[l,]
 pc.fdata<-sweep(pc.fdata,2,matrix(pc$mean$data,ncol=1),"+")
 basis.pc = pc$rotation[l, ]
 rownames(basis.pc$data) <- paste("PC", l, sep = "")
# pc.fdata<-inprod.fdata(x,pc$rotation)%*%pc$rotation$data
# pc.fdata<-inprod.fdata(fdataobj,basis.pc)%*%basis.pc$data

 basisobj<-pc
 fdnames<- colnames(pc$x[,l])
 if (is.null(basis)) {
   pc.fdata<-fdata(pc.fdata,tt,rtt,fdataobj$names)
   out <- list(fdataobj.pc=pc.fdata,basis = basis.pc, x = pc$x, mean = pc$mean,
   fdataobj.cen = pc$fdataobj.cen,fdataobj = pc$fdataobj,l = l,
   rn=rn,type = "pc")
   class(out) <- "fdata.comp"
   }
 else {
      fdobj<- Data2fd(argvals = tt, y = t(pc.fdata),basisobj = basis)
      out<-list()
      out$harmonics<-fdobj
      colnames(out$harmonics$coefs)<-rownames(fdataobj$data)
      out$values<-pc$newd^2
      out$scores<-pc$x[,l]
      rownames(out$scores)<-rownames(fdataobj$data)
      out$rn<-rn
      out$varprop<-out$values[l]/sum(out$values)
      out$meanfd<- Data2fd(argvals = tt, y = pc$mean$data[1,],basisobj = basis)
#      fdobj$type = "pca"
      class(out) <- "pca.fd"
      }
 return(out)
 # min.basis vs min.np
}


#######################
#######################
create.pls.basis<-function(fdataobj,y,l=1:5,rn=0){
     pls<-fdata2pls(fdataobj,y,ncomp=max(l))
     basis=pls$rotation[l,]
     rownames(basis$data)<-paste("PLS",l,sep="")
#     if (length(l)==1)   x=as.matrix(pls$x[,l],ncol=1)
#     else    x=pls$x[,l]
return(list("basis"=basis,"x"=pls$x,"mean"=pls$mean,"rn"=rn,
"fdataobj.cen"=pls$"fdataobj.cen","fdataobj"=pls$fdataobj,
"l"=l,"type"="pls","y"=y))
}
#######################
#######################
create.raw.fdata=function (fdataobj, l = 1:ncol(fdataobj))
{
    return(list(basis =fdataobj[,l] , type = "raw"))
}
#######################
#######################


