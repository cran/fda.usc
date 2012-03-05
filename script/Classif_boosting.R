################################################################################
# BOOSTING
#	INFERENCE PHONDAT

data(phoneme)
mlearn<-phoneme[["learn"]]
glearn<-phoneme[["classlearn"]]
mtest<-phoneme[["test"]]
gtest<-phoneme[["classtest"]]

ldata<-list("y"=glearn,"x"=mtest,"x2"=mlearn)
basis1=create.pc.basis(mlearn,1:3)
basis2=create.pc.basis(mtest,1:2)
lbasis<-list("x"=basis1,"x2"=basis2)

par.list2<-list()
par.list2[[1]]<-list("classif"="tree","fdataobj"=mlearn,nrep=2)
par.list2[[2]]<-list("classif"="tree","fdataobj"=mlearn,nrep=2)
out3<-classif.boosting(glearn,par.list=par.list2)
out3$maxprob
sum(out3$pred2!=glearn)/250
out3$pond
out3$fit
par.list2<-list()
par.list2[[1]]<-list("classif"="tree","fdataobj"=mlearn,nrep=4)
out2<-classif.boosting(glearn,par.list=par.list2)
out2$maxprob

par.list2<-list()
par.list2[[1]]<-list("classif"="tree","fdataobj"=mlearn,nrep=2)
par.list2[[2]]<-list("classif"="knn","fdataobj"=mlearn,nrep=2)
out4<-classif.boosting(glearn,par.list=par.list2)
out4$maxprob
out3$maxprob
sum(out3$pred2!=glearn)/250
sum(out4$pred2!=glearn)/250

par.list2<-list()
par.list2[[1]]<-list("classif"="glm","fdataobj"=mlearn,nrep=1)
par.list2[[2]]<-list("classif"="knn","fdataobj"=mlearn,nrep=1)
par.list2[[3]]<-list("classif"="tree","fdataobj"=mlearn,nrep=2)
out5<-classif.boosting(glearn,par.list=par.list2)
out4$maxprob
out5$maxprob
out3$maxprob
sum(out3$pred2!=glearn)/250
sum(out4$pred2!=glearn)/250
sum(out5$pred2!=glearn)/250
sum(out3$group.est!=glearn)/250
sum(out4$group.est!=glearn)/250
sum(out5$group.est!=glearn)/250


####################################
### PREDICCION #####################
# same data
par.list2<-list()
par.list2[[1]]<-list("classif"="tree","fdataobj"=mlearn,nrep=2)
out2<-classif.boosting(glearn,par.list=par.list2)

par.list2[[1]]<-list("classif"="tree","fdataobj"=mlearn,nrep=10)
out3<-classif.boosting(glearn,par.list=par.list2)
out2$maxprob
out3$maxprob

par.list3<-list(mlearn)
pred2<-predict.classif.boosting(out2,par.list3)
pred3<-predict.classif.boosting(out3,par.list3)
sum(pred2!=glearn)/250
sum(pred3!=glearn)/250

# new data
par.list3<-list(mtest)
pred2<-predict.classif.boosting(out2,par.list3)
pred3<-predict.classif.boosting(out3,par.list3)
sum(pred2!=glearn)/250
sum(pred3!=glearn)/250

#univariate
dataf<-data.frame(glearn)
dat=list("df"=dataf,"x"=mlearn)
a1<-classif.tree(glearn~x,data=dat)
newdat<-list("x"=mtest)
p1<-predict.classif(a1,type="class")
sum(p1!=gtest)/250
p1<-predict.classif(a1,newdat,type="class")
sum(p1!=gtest)/250



# new data
par.list2<-list()
par.list2[[1]]<-list("classif"="glm","fdataobj"=mlearn,nrep=2)
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred5!=glearn)/250

dataf<-data.frame(glearn,sample(250))
names(dataf)[2]<-"ind"
dat=list("df"=dataf,"mlearn"=mlearn)
a1<-classif.glm(glearn~mlearn, data = dat)
newdat<-list("mlearn"=mtest)
p1<-predict.classif(a1,newdat)
table(gtest,p1)
sum(p1==gtest)/250


#	ESTIMATION knn
out1=classif.knn(glearn,mlearn,knn=c(3,5,7))
pred1=predict.classif(out1,mtest)
sum(pred1 != gtest)/length(gtest)

# new data
par.list2<-list()
par.list2[[1]]<-list("classif"="knn","fdataobj"=mlearn,nrep=1)
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred55!=glearn)/250

par.list2<-list()
par.list2[[1]]<-list("classif"="knn","fdataobj"=mlearn,nrep=1)
par.list2[[2]]<-list("classif"="glm","fdataobj"=mlearn,nrep=1)
par.list2[[3]]<-list("classif"="tree","fdataobj"=mlearn,nrep=2)
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest,mtest,mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred55!=glearn)/250

par.list2<-list()
par.list2[[1]]<-list("classif"="kgam","fdataobj"=mlearn,nrep=2,
control=list(maxit=3))
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest,mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred55!=glearn)/250


par.list2<-list()
par.list2[[1]]<-list("classif"="gsam","fdataobj"=mlearn,nrep=1)
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest,mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred55!=glearn)/250

par.list2<-list()
par.list2[[1]]<-list("classif"="kgam","fdataobj"=mlearn,nrep=3,
control=list(maxit=3))
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest,mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred55!=glearn)/250


par.list2<-list()
par.list2[[1]]<-list("classif"="gsam","fdataobj"=mlearn,nrep=1)
par.list2[[2]]<-list("classif"="kgam","fdataobj"=mlearn,nrep=2,
control=list(maxit=3))
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest,mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred55!=glearn)/250

dataf<-data.frame(glearn)
dat=list("df"=dataf,"x"=mlearn)
a1<-classif.gsam(glearn~s(x),data=dat)
newdat<-list("x"=mtest)
p1<-predict.classif(a1,newdat)
sum(p1!=gtest)/250

a1<-classif.kgam(glearn~x,data=dat,control=list(maxit=3))
p1<-predict.classif(a1,newdat)
sum(p1!=gtest)/250


par.list2<-list()
par.list2[[1]]<-list("classif"="gsam","fdataobj"=mlearn,nrep=1)
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest,mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred55!=glearn)/250

#	ESTIMATION kernel
h=2^(0:5)
out2=classif.kernel(glearn,mlearn,h=h)
pred2=predict.classif(out2,mtest)
sum(pred2!= gtest)/length(gtest)
dataf<-data.frame(glearn)
dat=list("df"=dataf,"x"=mlearn)
a1<-classif.tree(glearn~x,data=dat)
newdat<-list("x"=mtest)
p1<-predict.classif(a1,newdat,type="class")
sum(p1==gtest)/250

par.list2<-list()
par.list2[[1]]<-list("classif"="knn","fdataobj"=mlearn,nrep=2)
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred5!=glearn)/250

par.list2<-list()
par.list2[[1]]<-list("classif"="kernel","fdataobj"=mlearn,nrep=2)
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred5!=glearn)/250

 par.list2<-list()
par.list2[[1]]<-list("classif"="tree","fdataobj"=mlearn,nrep=4)
out2<-classif.boosting(glearn,par.list=par.list2)
par.list3<-list(mtest)
pred55<-predict.classif.boosting(out2,par.list3)
sum(pred5!=glearn)/250




