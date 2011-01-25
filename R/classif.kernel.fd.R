classif.kernel.fd <- function(fdataobj,group,h=NULL,metric=metric.lp,...) {
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
data<-fdataobj[["data"]]
nn=ncol(data)
numgr=nrow(data)
tt =fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
if (numgr!=(length(group))) {stop("ERROR IN THE DATA DIMENSIONS") }
if (!is.factor(group)) {group=as.factor(group)}
if (numgr!=length(group)) {stop("ERROR IN THE DATA DIMENSIONS") }
C<-match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("fdataobj", "group", "h","metric"), names(mf), 0L)
group=as.factor(group)
numg=nlevels(as.factor(group))
group.pred=array(0,dim=c(numgr))
prob.groupV2=array(0,dim=c(numg))
pr=1;maxk=0;prob=0
mm=data
if (m[4]==0) mdist=metric(fdataobj,fdataobj,...) #            metric
else mdist=metric(data,data,...) #            metric
diag(mdist)=NA
Y=array(0,dim=c(numg,numgr)) 
for (g in 1:numg) {
y=as.integer(group==levels(group)[g])
Y[g,]=as.integer(group==levels(group)[g])
}
if (is.null(h)) {
     mdist2=mdist
     diag(mdist2)=NA
     h1=apply(mdist2,1,min,na.rm=TRUE)
     h.min=median(h1)
     h.max=max(h1)
     h.min=min(quantile(mdist2,probs=0.05,na.rm=TRUE),h.min)
     h.manx=max(quantile(mdist2,probs=0.25,na.rm=TRUE),h.max)
     h=seq(h.min,h.max,len=51)
  }
lenh=length(h)
group.est=array(0,dim=c(length(h),numgr))
pgrup=array(0,dim=c(numg,numgr,length(h)))
misclassification=array(1,dim=c(1,length(h)))
pb=txtProgressBar(min=0,max=lenh,style=3)
for (k in 1:lenh) {
   setTxtProgressBar(pb,k-0.5)
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
misclassification[1,k]=sum(group.est[k,]!=group)/numgr
if (pr>misclassification[1,k]) {
   pr=misclassification[1,k]
   iknn=k
   prob=1-pr
   group.pred=group.est[k,]
   prob.group2=t(pgrup[,,k])
   h.opt=h[k]
   D=mdist
   diag(D)=0  }
setTxtProgressBar(pb,k)}
close(pb)
misclassification[1,k]=sum(group.est[k,]!=group)/numgr
for (g in 1:numg) {
     g2=levels(group)[g]
     prob.groupV2[g]=sum((group.pred==g2)&(group.pred==group))/sum(group==g2)
     }
prob.groupV2=t(as.matrix(prob.groupV2))
colnames(misclassification)=round(h,2)
rownames(misclassification)="prob="
rownames(prob.groupV2)="prob="
colnames(prob.groupV2)=levels(group)
colnames(prob.group2)<-levels(group)
if (lenh>1) {
  if (h.opt==min(h))  cat(" Warning: h.opt is the minimum bandwidth vaulue provided, range(h)=",range(h),"\n")
  else  if (h.opt==max(h))  cat(" Warning: h.opt is the maximum bandwidth vaulue provided, range(hh)=",range(h),"\n")
}
output<-list(fdataobj=fdataobj,group=group,group.est=as.factor(group.pred),
max.prob=prob,h.opt=h.opt,D=D,prob.classification=prob.groupV2,
prob.group=t(prob.group2),misclassification=misclassification,h=h,
C=C,m=m)
class(output)="classif.fd"
return(output)
}

