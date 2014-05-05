################################################################################
# classif.DD: Fits Nonparametric Classification Procedure Based on DD–plot
# File created by Manuel Oviedo de la Fuente  using code from paper:
# Li, J., P.C., Cuesta-Albertos, J.A. and Liu, R.
# DD--Classifier: Nonparametric Classification Procedure Based on DD-plot.
# Journal of the American Statistical Association (2012), Vol. 107, 737--753.
################################################################################
classif.DD <- function(group,fdataobj,depth="FM",classif="glm",w,
par.classif=list(),par.depth=list(),
control=list(verbose=FALSE,draw=TRUE,col=NULL,alpha=.25)){
 C<-match.call()
 if (!is.factor(group)) group<-factor(group)
  lev<-levels(group)
 group<-factor(group,levels=lev[which(table(group)>0)])
 if (is.null(control$verbose))  control$verbose<-FALSE
 if (is.null(control$draw))  control$draw<-TRUE
 if (is.null(control$fine))  control$fine<-100
 if (is.null(control$alpha))  control$alpha<-0.25
#  if (is.null(control$bg))  control$bg<-"white"
# if (is.null(control$gray.scale))  control$gray.scale=FALSE
 classif0<-classif
 if (classif=="MaxD") {
 if (!is.null(par.classif$pol) & control$verbose) print("Maximum depth is done using polynomial argument, pol=1")
 par.classif$pol<-1
 classif<-"DD1" 
 }
# if (classif=="knn" | classif=="np" | classif=="grm")   control$fine<-min(50,control$fine)
 func.clas<-list()
 draw<-control$draw

 par.Df<-list("df"=data.frame(group))
# group<-ifelse(group==lev[1],0,1)
 ng<-length(lev)
 n<-length(group)#nrow(fdataobj)
 Df<-matrix(NA,ncol=ng,nrow=n)
 ind<-matrix(NA,nrow=n,ncol=ng)
 nvec<-table(group)
 p<-nvec[1]/n
 if ( is.fdata(fdataobj) | is.matrix(fdataobj) | is.data.frame(fdataobj)) {
  var.name<-deparse(substitute(fdataobj))
  ldata<-list(fdataobj)

  names(ldata)<-var.name
 }
 else {
  if (is.null(names(fdataobj))) names(fdataobj)<-paste("var",1:length(fdataobj),sep="")
  ldata<-fdataobj  
  var.name<-names(fdataobj)
 }
 lenlista<-length(ldata)
 lendepth<-length(depth)
 model<-FALSE
 par.ldata<-list()
 isfdata<-is.fdata(fdataobj) 
 if (is.null(names(par.depth))) {
   for (il in 1:lenlista)     par.ldata[[il]]<-par.depth     # esto no puede ser
 }
 else { 
  if (isfdata)  par.ldata[[1]]<-par.depth       
  else    par.ldata<-par.depth   }     
# print(par.ldata)
 lenl<-1
 ng2<-ng
 nam2<-c()
 nam<-NULL
 integrated<-FALSE
#  print(par.depth)
 if (missing(w)) {
   if (depth[1] %in% c("RPp","FMp","modep"))   {
     model=FALSE
################################################################################
# integrated version modal
# print("entra  integradora")  
   multi<-TRUE
   depth0<-depth
#   if (is.null(par.depth$scale)) par.depth$scale<-TRUE
   integrated<-TRUE  
    if (depth[1]=="modep"){
     hq<-numeric(ng)
     ismdist<-is.matrix(par.depth$metric)
     if (is.null(par.depth$par.metric$dscale)) par.depth$par.metric$dscale=mean
     if (ismdist)    mdist<-par.depth$metric
     fdataobj<-ldata[[1]]
     par.depth$lfdata<-ldata
     for (i in 1:ng) {
      ind[,i]<-group==lev[i]
      nam<-c(nam,paste("depth ",lev[i],sep=""))
      par.depth$lfdataref<-c.ldata(ldata,ind[,i]) #seleccionamos los de el grupo i
      if (ismdist)  par.depth$metric<-mdist[,ind[,i]]    
      oo<-do.call("depth.modep",par.depth)        #call a la funcion depth
      Df[,i]<-oo$dep
      hq[i]<-oo$hq      
     }     
     w<- attributes(oo$mdist)$method
     par.depth$h<-hq
     } 
   if (depth[1] %in% c("FMp","RPp")){ 
    par.depth$lfdata<-ldata
    fdataori<-ldata
    nam.depth<-paste("depth.",depth[1],sep="")

    for (i in 1:ng) {
      ind[,i]<-group==lev[i]
      fdataori<-c.ldata(ldata,ind[,i]) 
      nam<-c(nam,paste("depth ",lev[i],sep=""))
      par.depth$lfdataref<-fdataori 
      oo<-do.call(nam.depth,par.depth)
      Df[,i]<-oo$dep
      if (depth[1]=="RPp")  {par.depth$proj<-oo$proj     }
     }  
     w<-oo$dfunc
   }
   nam2<- paste(paste(var.name,collapse=".",sep=""),".",depth,".",lev,sep="")  
   gest<-factor(lev[apply(Df,1,which.max)],levels=lev) # Maximum depth  
   depthl<-depth
   colnames(Df)<-nam2 
    par.ldata<-par.depth      
 }
   else{
      lenl<-lenlista
      w<-rep(1/lenlista,len=lenlista)
      ng2<-ng*lenlista
      model<-TRUE
     }     

  }
if (!integrated){
 lenpesos<-length(w)
 if (sum(w)!=1) stop("Incorrect w argument, w must sum to 1")
 if (any(w<0))  stop("Incorrect w argument, w must be a positive")
 if (lenlista!=lenpesos) stop("Incorrect w argument")
 if (lendepth==1) depthl<-rep(depth,len=lenlista)
 else depthl<-depth
 depth0<-depth    
for (idat in 1:lenlista) {
 fdataobj<-ldata[[idat]]
 depth<-depthl[idat]
 par.depth<-par.ldata[[idat]]
 nc<-ncol(fdataobj)
 x<-array(NA,dim=c(n,nc,ng))
 Df<-matrix(NA,ncol=ng,nrow=n)
 ind<-matrix(NA,nrow=n,ncol=ng)
 isfdata<-is.fdata(fdataobj)
 mnames<-ls(pattern="^mdepth.*",envir=as.environment("package:fda.usc"),all.names=TRUE)
 fnames<-ls(pattern="^depth.*",envir=as.environment("package:fda.usc"),all.names=TRUE) 
 mnames2<-ls(pattern="^mdepth.*",envir=.GlobalEnv,all.names=TRUE) 
 fnames2<-ls(pattern="^depth.*",envir=.GlobalEnv,all.names=TRUE)  
 mnames<-c(mnames,mnames2)
 fnames<-c(fnames,fnames2) 
 depth.long<-paste("mdepth.",depth,sep="")  
 if (depth.long %in% mnames & !isfdata) {
  multi<-TRUE
  par.depth$x<-fdataobj
  }
 else    {
  depth.long<-paste("depth.",depth,sep="")
  if (depth.long %in% fnames & isfdata) { 
       par.depth[["fdataobj"]]<-fdataobj
       multi=FALSE         }
  else stop("Incorrect depth function or data class object")
 }     
     
 if (depth %in% c("RHS","RP","RPD","RT")){
 if (is.null(par.depth$proj)) {
  d <- nc#-1
  u <- matrix(runif(d*25,-1,1),25,d)
  norm <- sqrt(rowSums(u*u))
  arg <- u/norm
  if (!multi & isfdata)  par.depth$proj<-fdata(arg,fdataobj$argvals,fdataobj$rangeval)
  else   par.depth$proj<-arg
 }
 }
 ismdist<-is.matrix(par.depth$metric)
 if (ismdist) {
   mdist<-par.depth$metric
 }  
 dmode<-c(depth.long=="depth.mode" | depth.long=="mdepth.mode")
  if (dmode) hq<-numeric(ng) 
 for (i in 1:ng) {
   ind[,i]<-group==lev[i]
   nam<-c(nam,paste("depth ",lev[i],sep=""))#,paste("depth ",paste(lev[-i],collapse=",")))
   if (multi)  par.depth$xx<-fdataobj[ind[,i],]
   else   par.depth$fdataori<-fdataobj[ind[,i],]
   if (ismdist) {
     par.depth$metric<-mdist[,ind[,i]]   
     par.depth$metric2<-mdist[ind[,i],ind[,i]]
    }
     oo<-do.call(depth.long,par.depth)
     Df[,i]<-oo$dep
     if (dmode) hq[i]<-oo$hq    
  }
 if (dmode) par.depth$h<-hq
 for (i in 1:length(var.name)) var.name[i]<-unlist(strsplit(var.name[i], "[$]"))[[1]]
 for (i in 1:length(var.name)) var.name[i]<-unlist(strsplit(var.name[i], "[[]"))[[1]]     
 if (model) {
    if (idat==1) Df2<-Df
    else  Df2<-cbind(Df2,Df)
    par.ldata[[idat]]<-par.depth
    par.Df[[paste(".",idat,sep="")]]<-fdata(Df)
    nam2<-c(nam2, paste(var.name[idat],".",depthl[idat] ,".",lev,sep=""))
  }
  else{
    if (idat==1) {
    Df2<-w[idat]*Df
    }
    else  Df2<-Df2+w[idat]*Df
    nam2<- paste(nam2,paste(var.name[idat],depthl[idat],".",sep=""),sep="")
     par.ldata[[idat]]<-par.depth
  }
}
 if (!model) nam2<- paste(nam2,lev,sep="")
 Df<-Df2
 gest<-factor(lev[apply(Df,1,which.max)],levels=lev) # Maximum depth
 nvecs<-c(nvec,nvec)
 k<-1
colnames(Df)<-nam2
}
if (draw) {
  if (is.null(control$col))     col1<-c(4,2,3,1,5:14)[1:ng]
  else col1<-control$col
  if (ng2>2) {
 if (control$verbose) {
    if (ng>2)  warning("DD-plot for more than 2 levels not implemented yet")
    else       warning("DD-plot for more than 1 depth function not implemented yet")
    }
    pch0<-c(21,24,22,23,21,24,22,23,21,24,22,23)[1:ng]   
    pch1<-pch0[group]   
    }
   else {
    rg<-range(Df)
    dff<-diff(rg)*.03
     minsq<-max(0,rg[1]-dff)
     maxsq<-rg[2]+dff 
    sq<-seq(minsq,maxsq,len=control$fine)
    vec<-data.frame(expand.grid(sq,sq))
    names(vec)<-nam2        
     pch0<-c(21,3) 
     pch1<-pch0[group]
#     fill1<-c(4,2)
#     mycols <- adjustcolor(palette("default"), alpha.f = control$alpha)
#     opal <- palette(mycols)                
      fill1<- adjustcolor(col1, alpha.f = control$alpha)
 }                       
     col3<-col1
     col2<-col1[group]    
}
switch(classif,
 DD1={ 
  if (ng2>2) {stop("DD-plot for more than 2 levels not available")}
  if (is.null(par.classif$pol)) {
  b<-unique(Df[Df[,1]!=0&Df[,2]!=0,1]/Df[Df[,2]!=0&Df[,1]!=0,2])
  m <- length(b)
  mis <- rep(0,m)
  mis <- sapply(b,MCR0,Df,ind)
  b0 <- min(b[which.min(mis)])
 }
 else b0<-par.classif$pol
 group.est<-factor(as.numeric(b0*Df[,1]<Df[,2]))
 levels(group.est)=lev
 incorrect<-group.est!=group
 mis<-mean(incorrect) 
 par.classif$pol<-b0
 if (draw & ng2<=2) {
     rg<-range(Df)
     dff<-diff(rg)*.1
     bb2<-ifelse(b0*vec[,1]>vec[,2],0,1)
     }          
},
DD2={
 if (ng2>2) {stop("DD-plot for more than 2 levels not available")}
 if (is.null(par.classif$pol)) {
 if (is.null(par.classif$nmax))   nmax=5000   #0
   else   nmax<-par.classif$nmax
  nsample<-choose(n,2)
  combs1<-NULL
  if (nsample<nmax |  nmax==0 )  combs1 <- t(combn(n,2))
  else {
   for(i in 1:nmax)      combs1<-rbind(combs1,sample(n,2))
   if (control$verbose & nmax<nsample) warning("The number of polynomials considered is too large ",nsample,", the function search is a subsample of size, nmax=50000")
     }
    mcrs <- apply(combs1,1, quad.fit0, dep=Df,ind=ind)
    ind0 <- combs1[which.min(mcrs),] 
  A <- matrix(c( Df[,1][ind0[1]],( Df[,1][ind0[1]])^2, Df[,1][ind0[2]],( Df[,1][ind0[2]])^2),byrow=TRUE,2,2)
  ww <- c(Df[,2][ind0[1]],Df[,2][ind0[2]])
  a0.2 <-try(solve(A,ww),silent=TRUE)
  if (class(a0.2)=="try-error") {   #new
    sv<-svd(A)
    a0.2<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u))%*%ww)
     if (control$verbose)  warning("Inverse of sigma computed by SVD")
    }    
  group.est<-factor(ifelse(sapply(Df[,1],RR,a=a0.2)<Df[,2],lev[2],lev[1])   )
  incorrect<-group.est!=group
  mis<-mean(incorrect) 
  a2<-a0.2
  if (control$verbose) cat("pol=",a2," misclassification=",mis,"\n")
  if (is.null(par.classif$tt)) ind.tt<-100
  else ind.tt<-c(100,90, 20,80,140,70,160,60,180)
  for (ii in ind.tt) {
 a02 <- optim(a0.2,AMCR0,dep=Df,ind=ind,tt=ii)$par
 group.est<-factor(ifelse(sapply(Df[,1],RR,a=a02)<Df[,2],lev[2],lev[1])   ) 
 incorrect<-group.est!=group
 mis.new<-mean(incorrect)
 if (control$verbose) cat("pol=",a02," tt=",ii," misclassif=",mis.new,"\n")
 if (mis>=mis.new) {a2<-a02;mis<-mis.new;}
 }
 }
else a2<-par.classif$pol 
 group.est<-factor(ifelse(sapply(Df[,1],RR,a=a2)<Df[,2],lev[2],lev[1])   )
 incorrect<-group.est!=group
 mis<-mean(incorrect) 

 par.classif$pol<-a2
 if (draw & ng2<=2)  {
      bb2<-ifelse(sapply(vec[,1],RR,a=a2)>vec[,2],0,1)
     }
 },
 DD3={
 if (ng2>2) {stop("DD-plot for more than 2 levels not available")}
 if (is.null(par.classif$pol)) {
   if (is.null(par.classif$nmax))   nmax=500
   else   nmax<-par.classif$nmax
  nsample<-choose(n,3)
  combs1<-NULL
  if (nsample<nmax |  nmax==0)  combs1 <- t(combn(n,3))
  else {
   for(i in 1:nmax)      combs1<-rbind(combs1,sample(n,3))
   if (control$verbose & nmax<nsample) warning("The number of polynomials considered is too large ",nsample,", the function search is a subsample of size, nmax=50000")
     }   
 mcrs <- apply(combs1,1,cubic.fit0, dep=Df,ind=ind)
 ind0<- combs1[which.min(mcrs),]
 A <- matrix(c( Df[,1][ind0[1]],( Df[,1][ind0[1]])^2,( Df[,1][ind0[1]])^3,
           Df[,1][ind0[2]],( Df[,1][ind0[2]])^2,( Df[,1][ind0[2]])^3,
           Df[,1][ind0[3]],( Df[,1][ind0[3]])^2,( Df[,1][ind0[3]])^3),byrow=TRUE,3,3)
  ww <- c(Df[,2][ind0[1]],Df[,2][ind0[2]],Df[,2][ind0[3]])
  a0.3 <-try(solve(A,ww),silent=TRUE)
  if (class(a0.3)=="try-error") {   
     sv<-svd(A)
     a0.3<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u))%*%ww)
     if (control$verbose) warning("Inverse of sigma computed by SVD")
    }
 group.est<-factor(ifelse(sapply(Df[,1],RR,a=a0.3)>Df[,2],lev[1],lev[2]))   
  incorrect<-group.est!=group
 mis<-mean(incorrect) 

 a2<-a0.3
 if (control$verbose) cat("pol=",a2," misclassification=",mis,"\n")
 if (is.null(par.classif$tt)) ind.tt<-100
 else ind.tt<-c(100,90,120,80,140,70,160,60,180)
 for (ii in ind.tt) {
  a02 <- optim(a0.3,AMCR0,dep=Df,ind=ind,tt=ii)$par
   group.est<-factor(ifelse(sapply(Df[,1],RR,a=a02)>Df[,2],lev[1],lev[2])   )
  mis.new<-mean(group.est!=group)
  if (control$verbose) cat("pol=",a02," tt=",ii," misclassif=",mis.new,"\n")
  if (mis>=mis.new) {a2<-a02;mis<-mis.new}
 }
 }
 else a2<-par.classif$pol 
 group.est<-factor(ifelse(sapply(Df[,1],RR,a=a2)<Df[,2],lev[2],lev[1])   ) 
 mis<-mean(group.est!=group)
 par.classif$pol<-a2
  if (draw & ng2<=2)  {  
     bb2<-ifelse(sapply(vec[,1],RR,a=a2)>vec[,2],0,1)
     }        
 },
 lda={
   dat<-data.frame(as.numeric(group)-1,Df)
   par.classif$x<-Df
   par.classif$grouping<-group
   func.clas<-do.call(classif,par.classif)
   group.est<-predict(func.clas)$class
   mis<-mean(group.est!=group)
   if (draw & ng2<=2) {  
     bb2<-as.integer(predict(func.clas,vec)$class)
    }   
  },
 qda={
   dat<-data.frame(as.numeric(group)-1,Df)
   par.classif$x<-Df
   par.classif$grouping<-group    
   func.clas<-do.call(classif,par.classif)
   group.est<-predict(func.clas)$class
   mis<-mean(group.est!=group)
   if (draw & ng2<=2) {  
       bb2<-as.integer(predict(func.clas,vec)$class)
     }
 },
 knn={
  dat<-data.frame(as.numeric(group)-1,Df)
  names(dat)<-c("group1",nam2)
  par.classif$fdataobj<-Df
  par.classif$group<-group
  if (is.null(par.classif$knn))  par.classif$knn<-seq(3,19,by=2)
  func.clas<-do.call("classif.knn",par.classif)
  group.est<-func.clas$group.est
  mis<-mean(group.est!=group)
  if (draw & ng2<=2) {
       bb2<-as.numeric(predict(func.clas,vec))}
 },
 np={
  dat<-data.frame(as.numeric(group)-1,Df)
  names(dat)<-c("group1",nam2)
  par.classif$fdataobj<-Df
  par.classif$group<-group
  func.clas<-do.call("classif.kernel",par.classif)
  group.est<-func.clas$group.est
  mis<-mean(group.est!=group)
  if (draw & ng2<=2) {
     bb2<-as.numeric(predict.classif(func.clas,vec))}
   },
 tree={
   dat<-data.frame(group,Df)
   names(dat)<-c("group1",nam2)
   par.classif$formula<-formula(paste("group1~",names(dat)[-1]))
   par.classif$data<-dat
   func.clas<-do.call("rpart",par.classif)
   group.est<-predict(func.clas,dat,type="class")
   mis<-mean(group.est!=group)
   if (draw & ng2<=2) {
     bb2<-predict(func.clas,vec,type="class")
     bb2<-ifelse(bb2==lev[1],1,2)     
    }

 },
glm={
  dat<-data.frame(Df)
  names(dat)<-nam2
   par.classif$fdataobj<-dat
    par.classif$group<-group
  func.clas<-do.call("classif.glm2boost",par.classif)
  group.est<-factor(func.clas$group.est,levels=lev)
  mis<-mean(group.est!=group)
   if (draw & ng2<=2) {
    bb2<-predict.classif(func.clas,list("df"=vec),type="class")
    bb2<-ifelse(bb2==lev[1],1,2)
   }
 },
 gam={
  dat<-data.frame(Df)
  names(dat)<-nam2  
  nam.classif<-paste("classif.",classif,sep="")
  par.classif$fdataobj<-dat
  par.classif$group<-group                    
  if (is.null(par.classif$family)) par.classif$family<-binomial()
  func.clas<-do.call("classif.gsam2boost",par.classif)
  group.est<-func.clas$group.est
  incorrect<-group.est!=group
  mis<-mean(incorrect)
  if (draw & ng2<=2) {
   bb2<-predict.classif(func.clas,list("df"=vec),type="class")
   bb2<-ifelse(bb2==lev[1],1,2)
  }
  }
  )
 if (draw) {
    incorrect<-group.est!=group
    if (is.null(control$main)) {
    if (model) tit<-paste("DD-plot(",paste(depth0,collapse="-"),",",classif0,")",sep="")
    else    tit<-paste("DD-plot(",paste(depth0,w,collapse="-",sep=""),",",classif0,")",sep="")
    }
    else tit<-control$main   
    incorrect<-group!=group.est
    if (ng2>2) {    
       bg1<-col1[group]
       bg1[!incorrect] <-0
#       if (!control$gray.scale)               col2<-col1[group]       
       col2<-col1[group]       
       pairs(Df,col=col2,bg=bg1,main=tit,pch=pch1,diag.panel= function(...) {  
       legend("bottomleft",box.col=0,legend=lev,pch=pch0,col=col1,cex=max(.7,min(1,3/ncol(Df))),pt.bg=0,title="Points",title.adj=0.1,xjust=.1,yjust=.1)
       legend("bottomright",box.col=0,legend=c("good","bad"),pch=c(1,16),col=c(1,1),cex=min(1,3/ncol(Df)),pt.bg=c(0,1),title="Classify",title.adj=0.1,xjust=.1,yjust=.1)      
       } )} 
    else {
#     mycols <- adjustcolor(palette("default"), alpha.f = control$alpha)
#     opal <- palette(mycols)     
#       if (control$gray.scale) {
#          col2<-rep(1,len=15)[1:ng]
#          col1<-gray(c(.75,.5),alpha=control$alpha)[1:ng]
#          col3<-col2 
#          fill1<-col1
#      } 
     image(sq,sq,matrix(bb2,control$fine),xlab=nam[1],ylab=nam[2],main=tit,col=fill1,ylim=c(minsq,maxsq))    
 #    palette("default")
     points(Df,col=col2,lwd=1,bg=col2,pch=pch1)
 #    palette(mycols) 
     sinpuntos<-TRUE
     zona<-c("topright","bottomright","topleft","bottomleft")
     izona<-1                                                                                               
     if (is.null(control$box.col)) control$box.col="black"    #borde legenda
     if (is.null(control$bg)) control$bg="white"
     while (sinpuntos) {
      le<- legend(zona[izona],title="Data points",legend=lev,col=col3,pt.bg=col3,border=1,horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0.1,yjust=0,plot=FALSE)            #
      le2<-legend(zona[izona],title="Class Rule ",fill=fill1,legend=lev,border=1,pt.cex=2,pt.bg=col2,horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0.1,yjust=0,plot=FALSE)        
      sinpuntos<-switch(izona,
           "1"={
             xleg<-le$rect$left
             xleg2<-le2$rect$left             
             yleg2<-le$rect$top-le$rect$h -le2$rect$h
             yleg<-le$rect$top-le$rect$h
             sinpuntos= any(Df[,1]>xleg & Df[,2]>yleg2)},
           "2"={
             xleg<-le$rect$left
             xleg2<-le2$rect$left             
             yleg2<-le$rect$top-le2$rect$h  
             yleg<-le$rect$top                      
             sinpuntos= any(Df[,1]>xleg & Df[,2]<yleg)},
           "3"={
             xleg<-le$rect$left
             xleg2<-le2$rect$left             
             yleg2<-le$rect$top-le$rect$h -le2$rect$h
             yleg<-le$rect$top-le$rect$h
             sinpuntos= any(Df[,1]<(xleg+le$rect$w) & Df[,2]>yleg2) },
           "4"={ 
             xleg<-le$rect$left
             xleg2<-le$rect$left             
             yleg2<-le$rect$top-le2$rect$h
             yleg<-le$rect$top
             sinpuntos= any(Df[,1]<(xleg+le$rect$w) & Df[,2]<yleg2)})              
             
      if (!sinpuntos | izona==4) {
           if (!sinpuntos & izona==4 ) {
             control$box.col=0   #borde legenda
             control$bg=0
        le<- legend(zona[1],title="Data points",legend=lev,col=col3,pt.bg=col3,border=1,horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0.1,yjust=0,plot=FALSE)            #
        le2<-legend(zona[1],title="Class Rule ",fill=fill1,legend=lev,border=1,pt.cex=2,pt.bg=col1,horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0.1,yjust=0,plot=FALSE)         
             xleg<-le$rect$left
             xleg2<-le2$rect$left             
             yleg2<-le$rect$top-le$rect$h -le2$rect$h
             yleg<-le$rect$top-le$rect$h                     
             }     
        sinpuntos<-FALSE
        }
      else   izona<-izona+1
     }        

#xjust=0,yjust=0,      #pt.cex=2,

     legend(xleg2,yleg2, title="Class Rule ",fill=fill1,legend=lev,border=1, pt.bg=col2,
     horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0,yjust=0,box.col=control$box.col,bg=control$bg)             
#     palette("default")
     legend(xleg,yleg,title="Data points",legend=lev,pch=pch0,col=col3,pt.bg=col3,
     border=1,horiz=FALSE,cex=0.8,title.adj=0.5,xjust=0,yjust=0,box.col=control$box.col,bg=control$bg)
    }
    }
prob.classification<-diag(table(group,group.est))/table(group)
output<-list("group.est"=group.est,"misclassification"=mis,
"prob.classification"=prob.classification,"dep"=Df,"depth"=depthl,
 "par.depth"=par.ldata,"classif"=classif,"par.classif"=par.classif,"w"=w,
 "group"=group,"fdataobj"=fdataobj,C=C,"fit"=func.clas,"model"=model,multi=multi,ldata=ldata)
class(output)=c("classif.DD","classif")
return(output)
}

###########################################################################################################################
############################################################################################################################
predict.classif.DD<-function(object,new.fdataobj=NULL,type="predictive",...){
 if (is.null(new.fdataobj)) return(object$group.est)
 depth<-object$depth
 par.depth<-object$par.depth
 multi<-object$multi                           
if (is.fdata(new.fdataobj)| is.matrix(new.fdataobj) | is.data.frame(new.fdataobj)) {  ldata<-list(new.fdataobj)}
else {   ldata<-new.fdataobj   }
lenlista<-length(ldata)
lendepth<-length(depth)
par.ldata<-object$par.depth
w<-object$w
lenpesos<-length(w)
depth<-depthl<-object$depth
ng2<-ncol(object$dep)
par.Df<-list()  
group<-object$group
nn<-length(group)
par.depth<-object$par.depth  
lev<-levels(group)
ng<-length(lev)
ind<-matrix(NA,nrow=nn,ncol=ng)
integrated<-FALSE   
 if (missing(w)) {   
    model<-TRUE
   lenl<-lenlista
   w<-rep(1/lenlista,len=lenlista)
   ng2<-ng*lenlista
 }
 else {
   if (depth[1]=="modep") {
    integrated<-TRUE
    n<-nrow(ldata[[1]])
    Df<-matrix(NA,ncol=ng,nrow=n)
    if (is.null(par.depth$par.metric$dscale)) par.depth$par.metric$dscale=mean
    hq<-par.depth$h
    par.depth$lfdata<-new.fdataobj
    for (i in 1:ng) {
      ind[,i]<-group==lev[i]
      par.depth$lfdataref<-c.ldata(object$par.depth$lfdata,ind[,i])
      par.depth$h<-hq[i]
      aa<-do.call("depth.modep",par.depth)$dep        #call a la funcion depth
      Df[,i]<-aa
     }    
   }
   if (depth[1]=="FMp" | depth[1]=="RPp"){ 
    integrated<-TRUE
    n<-nrow(new.fdataobj[[1]])
    Df<-matrix(NA,ncol=ng,nrow=n)       
    par.depth$lfdata<-new.fdataobj
    nam.dep<-paste("depth.",depth[1],sep="")
    for (i in 1:ng) {
      ind[,i]<-group==lev[i]
      par.depth$lfdataref<-c.ldata(object$par.depth$lfdata,ind[,i])
        Df[,i]<-do.call(nam.dep,par.depth)$dep        
     }    
   }  
}     
if (!integrated){
for (idat in 1:lenlista) {
 depth<-depthl[idat]
 par.depth<-par.ldata[[idat]]
 fdataobj<-par.depth$fdataobj
 new.fdataobj<-ldata[[idat]]
 nc<-ncol(fdataobj)
 if (multi) {
  depth.long<-paste("mdepth.",depth,sep="")
  if (is.vector(new.fdataobj)) new.fdataobj<-par.depth$x<-matrix(new.fdataobj,nrow=1)
  else par.depth$x<-new.fdataobj
  fdataobj<-object$par.depth[[idat]]$x
 }
 else    {
    depth.long<-paste("depth.",depth,sep="")
    if (is.vector(new.fdataobj[["data"]])) new.fdataobj[["data"]]<-matrix(new.fdataobj,nrow=1)
    par.depth[["fdataobj"]]<-new.fdataobj
    fdataobj<-object$par.depth[[idat]]$fdataobj
         }    
 n<-nrow(new.fdataobj)
 nc<-ncol(new.fdataobj)
 nvec<-table(group)
 p<-nvec[1]/n    
 Df<-matrix(NA,ncol=ng,nrow=n) 
 ismdist<-is.matrix(par.depth$metric)
 if (ismdist) {
   mdist<-par.depth$metric
   par.depth$metric<-attr(mdist, "call")
 }
 dmode<-c(depth.long=="depth.mode" | depth.long=="mdepth.mode")
 if (dmode)  hq<-par.depth$h
 for (i in 1:ng) {
    ind[,i]<-group==lev[i]
    if (multi)  par.depth$xx<-fdataobj[ind[,i],]
    else   par.depth$fdataori<-fdataobj[ind[,i],]
    if (dmode)      par.depth$h<-hq[i]
     Df[,i]<-do.call(depth.long,par.depth)$dep
 } 
if (object$model) {
    if (idat==1) Df2<-Df
    else  Df2<-cbind(Df2,Df)
    par.ldata[[idat]]<-par.depth
    par.Df[[paste("dep",idat,sep="")]]<-fdata(Df)
  }
  else{
    if (idat==1) Df2<-w[idat]*Df
    else  Df2<-Df2+w[idat]*Df
    par.ldata[[idat]]<-par.depth
  }
 } 
 Df<-Df2
 }
colnames(Df)<-colnames(object$dep)
group.est<-switch(object$classif,
# MD={gest<-apply(Df,1,which.max)},
 DD1={
   group.est<- factor(ifelse(object$par.classif$pol*Df[,1]>Df[,2],lev[1],lev[2]),levels=lev)},   
 DD2={
 group.est<- factor(ifelse(sapply(Df[,1],RR,a=object$par.classif$pol)>Df[,2],lev[1],lev[2]),levels=lev)
 },
 DD3={
  group.est<- factor(ifelse(sapply(Df[,1],RR,a=object$par.classif$pol)>Df[,2],lev[1],lev[2]),levels=lev)
  },
  
 lda={predict(object$fit,Df)$class},
 qda={predict(object$fit,Df)$class},
 glm={
   dat<-data.frame(Df)
    group.est<-predict.classif(object$fit,list("df"=dat),type = "class")
 },
 gam={
  dat<-data.frame(Df)
  group.est<-predict.classif(object$fit,list("df"=dat),type = "class")
  },
 tree={
  dat<-data.frame(Df)
  names(dat)<-colnames(object$dep)
  group.est<-predict(object$fit,dat,type = "class")
 },
 knn={
  dat<-data.frame(Df)
     group.est<-predict.classif(object$fit,dat,type = "class")
 },
 np={
  dat<-data.frame(Df)
  group.est<-predict(object$fit,dat)
 },
  grm={
  dat<-data.frame(Df)
  group.est<-predict.classif(object$fit,par.Df,type = "class")
  }
 )
 if (type=="dep") return(list("group.pred"=group.est,"dep"=Df))
 else   return(group.est)
}
