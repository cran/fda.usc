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
create.pc.basis<-function(fdataobj,l=1:5){
     pc<-fdata2pc(fdataobj,ncomp=max(l))
     basis=pc$rotation[l,]
     rownames(basis$data)<-paste("PC",l,sep="")
     if (length(l)==1)   x=as.matrix(pc$x[,l],ncol=1)
     else    x=pc$x[,l]
return(list("basis"=basis,"x"=x,"mean"=pc$mean,"type"="pc"))
}
#######################
#######################
create.raw.fdata=function (fdataobj, l = 1:ncol(fdataobj))
{
    return(list(basis =fdataobj[,l] , type = "raw"))
}
#######################
#######################
  #######################
create.pls.basis<-function(fdataobj,y,l=1:5){
     pls<-fdata2pls(fdataobj,y,ncomp=max(l))
     basis=pls$rotation[l,]
     rownames(basis$data)<-paste("PLS",l,sep="")
     if (length(l)==1)   x=as.matrix(pls$x[,l],ncol=1)
     else    x=pls$x[,l]
return(list("basis"=basis,"x"=x,"mean"=pls$mean,"type"="pls","y"=y))
}

