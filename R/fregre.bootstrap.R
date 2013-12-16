fregre.bootstrap<-function(model,nb=500,wild=TRUE,type.wild="golden",newX=NULL,smo=0.1,smoX=0.05,alpha=0.95,
kmax.fix=FALSE,draw=TRUE,...){
fdataobj=model$fdataobj
nas<-apply(fdataobj$data,1,count.na)
if (any(nas))  {
   fdataobj$data<-fdataobj$data[!nas,]
   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
   }
dat<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
nam<-fdataobj[["names"]]
resi=model$residuals
if (model$call[[1]]=="fregre.pc") {
     beta.est=model$beta.est$data
     pc=1
     }
else if (model$call[[1]]=="fregre.pls") {
          beta.est=model$beta.est
          pc=2
          }
else if (model$call[[1]]=="fregre.basis" || model$call[[1]]=="fregre.basis.cv") {
          beta.est=model$beta.est
          beta.est=eval.fd(tt,beta.est)
          pc=3
          }
else stop("No fregre.pc, fregre.basis or fregre.basis.cv object in model argument")
a.est=model$coefficients[1]
sr2=model$sr2
n <- nrow(fdataobj)
J <- ncol(fdataobj)
cb.num <- round(alpha * nb)
betas.boot <- array(NA,dim=c(nb,J))
betas.boot2<-model$beta.est
norm.boot <- array(NA,dim=c(nb,1))
ncoefs<-100

pb=txtProgressBar(min=1,max=nb,width=50,style=3)
knn.fix=NULL
y.mue2<-array(NA,dim=c(nb,nrow(dat)))
ypred<-array(NA,dim=c(nb,nrow(newX)))
  if (!is.logical(kmax.fix)) {
  criteria=kmax.fix
  kmax.fix=TRUE    } 
  else   {   if (pc<3)      criteria="SIC"
             if (pc==3)     criteria=GCV.S}
  if (pc==3) ncoefs<-nrow(model$beta.est$coefs)
  if (kmax.fix) coefs.boot <- array(NA,dim=c(nb,ncoefs))
  else coefs.boot<-list()
if (!wild){           
for (i in 1:nb){
   setTxtProgressBar(pb,i-0.5)
   muee <- sample(1:n,n,replace=TRUE)
   mueX <- sample(1:n,n,replace=TRUE)
   residuals.mue <- resi[muee] + rnorm(n,0,sqrt(smo * sr2))  
   b1<- fdata(mvrnorm(n,rep(0,J),smoX * var(dat)),argvals(fdataobj),rtt)
   b0<-fdataobj[mueX,]
   fdata.mue <-b0+b1
   if (pc==1)   {
      y.mue<-predict.fregre.fd(model,fdata.mue)  + residuals.mue
      if (kmax.fix)    funcregpc.mue <- fregre.pc(fdata.mue,y.mue,l=model$l,lambda=model$lambda,weights=model$weights,...)    
       else     {
               fpc <- fregre.pc.cv(fdata.mue,y.mue,max(model$l,8),lambda=model$lambda,criteria=criteria,weights=model$weights,...)
               knn.fix[[i]]<-fpc$pc.opt
               funcregpc.mue<-fpc$fregre.pc
                        }
       betas.boot[i,] <- funcregpc.mue$beta.est$data
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i] <- norm.fdata(bb)
       if (!is.null(newX))   {
       ypred[i,]<-predict(funcregpc.mue,newX)
       y.mue2[i,]<-y.mue
                               }
             }
else  if (pc==2)  {
      y.mue<-predict.fregre.fd(model,fdata.mue)  + residuals.mue

      if (kmax.fix)    {funcregpc.mue <- fregre.pls(fdata.mue,y.mue,model$l,...)}
       else     {
                        fpc <- fregre.pls.cv(fdata.mue,y.mue,max(model$l,8),criteria=criteria,...)
                        knn.fix[[i]]<-fpc$pls.opt
                        funcregpc.mue<-fpc$fregre.pls
                        }
       betas.boot[i,] <- funcregpc.mue$beta.est$data
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i] <- norm.fdata(bb)
       if (!is.null(newX))   {
       ypred[i,]<-predict(funcregpc.mue,newX)
       y.mue2[i,]<-y.mue
                               }
             }
  else  {
       bett<-fdata(t(beta.est),tt,rtt)
        y.mue<-predict.fregre.fd(model,fdata.mue)  + residuals.mue
       if (kmax.fix) funcregpc.mue <-fregre.basis(fdata.mue,y.mue,model$basis.x.opt,
        model$basis.b.opt,Lfdobj=model$Lfdobj,weights=model$weights,...)
       else {
           funcregpc.mue <-fregre.basis.cv(fdata.mue,y.mue,type.CV=criteria,Lfdobj=model$Lfdobj,weights=model$weights,...)                   
              knn.fix[[i]]<-c(funcregpc.mue$basis.x.opt$nbasis,funcregpc.mue$basis.b.opt$nbasis)
                        }
       betas.boot[i,] <- eval.fd(tt,funcregpc.mue$beta.est)
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i]<-  norm.fd(bb)
        if (kmax.fix)  coefs.boot[i,]<-funcregpc.mue$beta.est$coefs[,1]
        else         coefs.boot[[i]]<-funcregpc.mue$beta.est$coefs[,1]
       if (!is.null(newX))   {
       ypred[i,]<-predict(funcregpc.mue,newX)
       y.mue2[i,]<-y.mue      
                               }
       }      
   setTxtProgressBar(pb,i)             }    
close(pb)
}
else {
pred<-model$fitted.values
fdata.mue<-fdataobj
for (i in 1:nb){
   setTxtProgressBar(pb,i-0.5)
   muee <- sample(1:n,n,replace=TRUE)
   residuals.mue <- rwild(resi[muee],type.wild)
   fdata.mue <- fdataobj[muee] 
   if (pc==1)   {
      y.mue<-pred  + residuals.mue
      if (kmax.fix)    funcregpc.mue <- fregre.pc(fdata.mue,y.mue,l=model$l,lambda=model$lambda,weights=model$weights,...)
       else     {
               fpc <- fregre.pc.cv(fdata.mue,y.mue,max(model$l,8),lambda=model$lambda,P=model$P,criteria=criteria,weights=model$weights,...)
               knn.fix[[i]]<-fpc$pc.opt
               funcregpc.mue<-fpc$fregre.pc
                }
       betas.boot[i,] <- funcregpc.mue$beta.est$data
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i] <- norm.fdata(bb)
       if (!is.null(newX))   {
         ypred[i,]<-predict(funcregpc.mue,newX)
         y.mue2[i,]<-y.mue         }
       }
   else  if (pc==2)  {
      y.mue<-pred + residuals.mue
      if (kmax.fix)    {funcregpc.mue <- fregre.pls(fdata.mue,y.mue,model$l,...)}
       else     {
                        fpc <- fregre.pls.cv(fdata.mue,y.mue,max(model$l,8),criteria=criteria,...)
                        knn.fix[[i]]<-fpc$pls.opt
                        funcregpc.mue<-fpc$fregre.pls
                        }
       betas.boot[i,] <- funcregpc.mue$beta.est$data
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i] <- norm.fdata(bb)
       if (!is.null(newX))   {
         ypred[i,]<-predict(funcregpc.mue,newX)
         y.mue2[i,]<-y.mue   }
      }
    else  {
       bett<-fdata(t(beta.est),tt,rtt)
       y.mue<-pred + residuals.mue
       if (kmax.fix) funcregpc.mue <-fregre.basis(fdata.mue,y.mue,model$basis.x.opt,
        model$basis.b.opt,Lfdobj=model$Lfdobj,weights=model$weights,...)
       else {
           funcregpc.mue <-fregre.basis.cv(fdata.mue,y.mue,type.CV=criteria,Lfdobj=model$Lfdobj,weights=model$weights,...)                   
           knn.fix[[i]]<-c(funcregpc.mue$basis.x.opt$nbasis,funcregpc.mue$basis.b.opt$nbasis)
            }
       betas.boot[i,] <- eval.fd(tt,funcregpc.mue$beta.est)
       bb<-model$beta.est-funcregpc.mue$beta.est
       norm.boot[i]<-  norm.fd(bb)
       if (kmax.fix)  coefs.boot[i,]<-funcregpc.mue$beta.est$coefs[,1]
       else         coefs.boot[[i]]<-funcregpc.mue$beta.est$coefs[,1]
       if (!is.null(newX))   {
        ypred[i,]<-predict(funcregpc.mue,newX)
       y.mue2[i,]<-y.mue      
       }
       }      
   setTxtProgressBar(pb,i)             }    
close(pb)
}             
betas.boot<- fdata(betas.boot,tt,rtt,nam)
betas.boot$names$main<-"beta.est bootstrap"
if (draw) {
  out<-norm.boot>quantile(norm.boot,alpha)
  plot(betas.boot[-out],col="grey")
  lines(model$beta.est,col=4)
  lines(betas.boot[out],col=2,lty=2)
if (!is.null(newX))   {
   dev.new()
   IC<-apply(ypred,2,quantile,c((1-alpha)/2,alpha+(1-alpha)/2))
   yp<-predict(model,newX)
   if (nrow(newX)>1) matplot(t(rbind(IC,yp)),type="l",col=c(2,2,4),lty=c(2,2,1),
   xlab="Id curves",ylab="predicted value",main="y predicted and CI")
   else   plot(rep(1,3),c(IC,yp),type="p",col=c(2,2,4),lty=c(2,2,1),
   xlab=rownames(newX$data),xlim=c(0.9,1.1),ylab="predicted value",main="y predicted and CI")
}
}
return(list("betas.boot"=betas.boot,"norm.boot"=norm.boot,"coefs.boot"=coefs.boot,
"knn.fix"=knn.fix,"ypred"=ypred,"y.boot"=y.mue2))
}





