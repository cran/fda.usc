#######################
#######################
classif.knn<-function(group,fdataobj,w=NULL,
knn=seq(3,floor(min(table(group))/3),by=2),metric=metric.lp,...) {
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
if (numgr!=(length(group)))  {stop("ERROR IN THE DATA DIMENSIONS") }
if (min(table(group))<max(knn)) {
 cat(paste("Warning: the  class",which.min(table(group))," has only",
 min(table(group)),"elements or components"),"\n")}
C<-match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c( "group","fdataobj","w","knn","metric"), names(mf), 0L)
numg=nlevels(as.factor(group))
misclassification=array(1,dim=c(1,length(knn)))
group.pred=array(0,dim=c(numgr))
group.est=array(0,dim=c(length(knn),numgr))
prob.groupV2=array(0,dim=c(numg))
prob.group=array(0,dim=c(numgr,numg,length(knn)))
pr=1;
mm=data
mdist=array(0,dim=c(numgr,numgr))
mat=array(0,dim=c(numgr,numg,length(knn)))
if (m[5]==0) mdist=metric(fdataobj,fdataobj,...)
else mdist=metric(data,data,...)
pb=txtProgressBar(min=0,max=numgr,style=3)
for (i in 1:numgr) {
   setTxtProgressBar(pb,i-0.5)
   for (k in 1:length(knn)) {
     kk=knn[k]+1
     vmax=vvmax=rep(0,numgr)
     mdist[i,i]=NA
     vec=quantile(mdist[i,],prob=(kk/(numgr-1)),type=4,na.rm=TRUE)
     l=which(mdist[i,]<=vec)
     l2=which.min(mdist[i,])
     for (j in 1:knn[k]) {
         mat[i,group[l[j]],k]=mat[i,group[l[j]],k]+1
     }
     prob.group[i,,k]=mat[i,,k]/knn[k]
     max2=which(max(mat[i,,k])==mat[i,,k])
     if (length(max2)>1) {vvmax[i]=as.character(group[l2])}
     else {vvmax[i]=levels(group)[max2]}
     group.est[k,i]=vvmax[i]
   }
   setTxtProgressBar(pb,i)
   }
close(pb)
for (k in 1:length(knn)) {
   misclassification[1,k]=sum(group.est[k,]!=group)/numgr
   if (pr>misclassification[1,k]) {
    pr=misclassification[1,k]
    knn.opt=knn[k]
    max.prob=1-pr
    group.pred=group.est[k,]
    prob.group2=prob.group[,,k]
    D<-mdist
    diag(D)=0
   }
   for (g in 1:numg) {
    g2=levels(as.factor(group))[g]
    prob.groupV2[g]=sum((group.pred==g2)&(group.pred==group))/sum(group==g2)
   }
  }
  prob.groupV2=t(as.matrix(prob.groupV2))
 	rownames(prob.groupV2)="prob="
	colnames(prob.groupV2)=levels(group)
  colnames(prob.group2)<-levels(group)
  colnames(misclassification)=knn
  rownames(misclassification)="prob="
output<-list(fdataobj=fdataobj,group=group,group.est=as.factor(group.pred),
max.prob=max.prob,knn.opt=knn.opt,D=D,prob.classification=prob.groupV2,
prob.group=prob.group2,misclassification=misclassification,knn=knn,
C=C,m=m)
class(output)="classif"
return(output)
}


