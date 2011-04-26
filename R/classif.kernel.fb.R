classif.kernel.fb<-function(fdataobj,group,h=NULL,metric=metric.lp,
type.basis="bspline",par.basis=list(rangeval=NULL,nbasis=NULL),...)
{
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
if (is.null(names(group))) names(group)<-1:length(group)
nas<-apply(fdataobj$data,1,count.na)
nas.g<-is.na(group)
if (any(nas) & !any(nas.g)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
   group<-group[bb]
   }
else {
if (!any(nas) & any(nas.g)) {
   cat("Warning: ",sum(nas.g)," values of group with NA are omited \n")
   bb<-!nas.g
   fdataobj$data<-fdataobj$data[bb,]
   group<-group[bb]
   }
else {
if (any(nas) & any(nas.g))  {
   bb<-!nas & !nas.g
   cat("Warning: ",sum(!bb)," curves  and values of group with NA are omited \n")
   fdataobj$data<-fdataobj$data[bb,]
   group<-group[bb]
   }
}}
data<-fdataobj$data
nn=ncol(data)
numgr=nrow(data)
tt =fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
if (numgr!=(length(group))) { stop("ERROR IN THE DATA DIMENSIONS") }
if (!is.factor(group)) {group=as.factor(group)}
if (numgr!=(length(group)))   {stop("ERROR IN THE DATA DIMENSIONS") }
group=as.factor(group)
numg=nlevels(as.factor(group))
group.pred=array(0,dim=c(numgr))
prob.groupV2=array(0,dim=c(numg))
pr=1;maxk=0;prob=0
mm=data
if (is.null(h)) lenh=21
else lenh=length(h)
if (is.null(par.basis$nbasis)) nbas=floor(seq(from=dim(fdataobj)[1]/4,to=dim(fdataobj)[1]/8,len=10))
else nbas=par.basis$nbasis
if (is.null(par.basis$rangeval)) par.basis$rangeval=rtt
lenb=length(nbas)
pgrup=array(0,dim=c(numg,numgr,lenh))
misclassification=array(1,dim=c(lenb,lenh))
group.est=array(NA,dim=c(lenh,numgr))
pb=txtProgressBar(min=0,max=lenb,style=3)
for (bb in 1:lenb) {
   setTxtProgressBar(pb,bb-0.5)
    if (type.basis=="fourier") {
       if ((nbas[bb]%%2)==0) {nbas[bb]=nbas[bb]+1}}
   par.basis$nbasis<-nbas[bb]
 C <- match.call()                            #
 mf <- match.call(expand.dots = FALSE)        #
 m<-match(c("fdataobj", "group","h","metric","type.basis","par.basis"),names(mf),0L)
 imetric <- m[6]
 if (imetric == 0)  a1 <- create.bspline.basis
 else       a1 <- paste('create.',type.basis,'.basis',sep="")
   basis.obj<-do.call(a1,par.basis)
   mdist=array(0,dim=c(numgr,numgr))
   mm.fd=Data2fd(y=t(mm),argvals=tt,basisobj=basis.obj)
   x1<-t(mm.fd$coefs)
   fdata(x1,tt)

#if (m[4]==0) mdist=metric(fdataobj,fdataobj,...) #            metric
#if (m[4]==0) #mdist=metric(fdata(x1,tt),fdata(x1,tt),...) #            metric
# else mdist=metric(x1,x1,...) #            metric
mdist=metric(x1,x1,...)
if (is.null(h)) {
     mdist2=mdist
     diag(mdist2)=NA
     h1=apply(mdist2,1,min,na.rm=TRUE)
     h.min=median(h1)
     h.max=max(h1)
     h.min=min(quantile(mdist2,probs=0.025,na.rm=TRUE),h.min)
     h.manx=max(quantile(mdist2,probs=0.25,na.rm=TRUE),h.max)
     h=seq(h.min,h.max,len=lenh)
  }
   diag(mdist)=NA
   for (k in 1:lenh) {
        vmax=array(0,dim=c(numgr));vvmax=array(0,dim=c(numgr))
        mat=array(0,dim=c(numgr,numg))
        hmdist=mdist/h[k]
        for (j in 1:numg) {
            y=as.integer(group==levels(group)[j])
            pgrup[j,,k]<-apply(hmdist,1,rkernel,y=y)

            }
        for (i in 1:numgr) {
            l=which.max(pgrup[,i,k])
            group.est[k,i]=levels(group)[l]
        }
        misclassification[bb,k]=sum(group.est[k,]!=group)/numgr
        if (pr>misclassification[bb,k]) {
           pr=misclassification[bb,k]
           h.opt=h[k]
           iknn=k
           prob=1-pr
           group.pred=group.est[k,]
           prob.group2=pgrup[,,k]
           basis.opt=nbas[bb]
           basis.obj.opt=basis.obj
           coefs=mm.fd$coefs
           D=mdist
           diag(D)=0
           }
   }
   setTxtProgressBar(pb,bb)
   }
close(pb)
misclassification[bb,k]=sum(group.est[k,]!=group)/numgr
for (g in 1:numg) {
      g2=levels(group)[g]
      prob.groupV2[g]=sum((group.pred==g2)&(group.pred==group))/sum(group==g2)
   }
prob.groupV2=t(as.matrix(prob.groupV2))
dimnames(misclassification)[1]=list(nbas)
dimnames(misclassification)[2]=list(round(h,2))
colnames(prob.groupV2)=levels(group)
rownames(prob.groupV2)="prob="
dimnames(prob.group2)[1]=list(levels(group))
if (lenh>1) {
  if (h.opt==min(h))  cat(" Warning: h.opt is the minimum bandwidth vaulue provided, range(h)=",range(h),"\n")
  else  if (h.opt==max(h))  cat(" Warning: h.opt is the maximum bandwidth vaulue provided, range(hh)=",range(h),"\n")
}
output<-list(fdataobj=fdataobj,group=group,group.est=as.factor(group.pred),
max.prob=prob,h.opt=h.opt,D=D,prob.classification=prob.groupV2,
prob.group=t(prob.group2),basis.opt=basis.opt,basis.obj.opt=basis.obj.opt,
coefs=coefs,misclassification=misclassification,num.basis=nbas,
h=h,C=C,m=m)
class(output)="classif.fd"
return(output)
}
