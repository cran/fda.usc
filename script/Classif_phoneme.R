################################################################################
# We present three non parametric methods for discriminating functional data
# 1.k-Nearest Neighbor method for functional data
# 2.Non parametric kernel method for functional data
# 3.Non parametric kernel method for basis representation of functional data
# For more details see:
# -Ferraty, F. and Vieu, P. (2006). Nonparametricc functional data analysis.
#  Springer Series in Statistics, New York.
# -Ferraty, F. and Vieu, P. (2006). NPFDA in practice.
#  Free access on line at http://www.lsp.ups-tlse.fr/staph/npfda/
################################################################################

################################################################################
# Read training data
################################################################################
data(phoneme)
# Read learn curves
mlearn<-phoneme[["learn"]]
# Read learn groups
glearn<-phoneme[["classlearn"]]

################################################################################
#	ESTIMATION (output of class "classif.fd")
# We use the training sample, which is  known for each data functional group
# You can use different metric  functions or types of distance metric
# by changing the parameters: p and w of metric.lp function
################################################################################
#1. the argument "knn" is vector of number of nearest neighbors considered
k<-seq(3,11,by=2)
out1=classif.knn(glearn,mlearn,knn=k)

#2. the argument "h" is vector of bandwidth considered
h=h.default(mlearn,prob=c(0.001,0.1),len=10)
out2=classif.kernel(glearn,mlearn,h=h)


################################################################################
# Summarizes information from classification methods of class "classif.fd", shows:
# -Probability of correct classification by group: prob.classification.
# -Confusion matrix between the theoretical groups and estimated groups.
# -Highest probability of correct classification max.prob.
################################################################################
#1. If the object is returned from the function: classif.knn.fd()
# -Vector of probability of correct classification by number of neighbors knn.
# -Optimal number of neighbors: knn.opt.
summary(out1)

#2. If the object is returned from the function: classif.kernel.fd
#   -Vector of probability of correct classification by banwidth h.
#  -Functional measure of closeness (optimal distance, h.opt).
summary(out2)

################################################################################
# Read test curves
mtest<-phoneme[["test"]]

#	PREDICTION  classmember from new curves
pred1=predict(out1,mtest)
names(pred1)
#estimated groups: pred1$group.pred
#for each curve shows the probability of each group: pred1$prob.group
pred2=predict(out2,mtest)

# Read theoretical groups
gtest<-phoneme[["classtest"]]

# Confusion matrix between the theoretical groups and estimated groups
t1<-table(pred1,gtest)
t2<-table(pred2,gtest)

# Missclassification  rate
ntest<-length(gtest)
MisclasPredict1 <- sum(pred1!= gtest)/ntest
MisclasPredict2 <- sum(pred2 != gtest)/ntest

c(MisclasPredict1,MisclasPredict2)

# Probability of correct classification by group: prob.classification.
prob.class<-cbind(diag(t1),diag(t2))/50
#rownames(prob.class)<-c("aa","ao","dcl","iy","sh")
colnames(prob.class)<-c("knn.fd","kernel.fd")
prob.class
matplot(prob.class,ylab="Probability Correct Classification",xlab="Phonemes")



# ejemplo Cuevas Bailo clasificacion 2 primeros grupos por CV
#Data set k-NNj1 k-NNj2    PLS h-modal RP(hM)    MWR
#Phoneme  0.7300 0.7800 0.7400  0.7300 0.7450 0.6950

data(phoneme)
ind<-151:250
mlearn<-phoneme[["learn"]][ind]
glearn<-phoneme[["classlearn"]][ind]
mtest<-phoneme[["test"]][ind]
gtest<-phoneme[["classtest"]][ind]

x<-c(mlearn,mtest)
y<-factor(c(glearn,glearn))
kk<-ng<-length(y)
pred<-matrix(NA,nrow=9,ncol=ng)
#	ESTIMATION knn
for (i in 1:ng) {
 cat(i)
 xx<<-x[-i,]
 yy<<-y[-i]
# out1=classif.knn(yy,xx)
# out2=classif.kernel(yy,xx)
# out3=classif.kgam2boost(yy,xx)
# out4=classif.gsam2boost(yy,xx)
# out5=classif.glm2boost(yy,xx)
# out6=classif.tree2boost(yy,xx)
 #kk[i]<-out1$knn.opt
 newx<-x[i]
# pred[1,i]=predict.classif(out1,newx)
# pred[2,i]=predict.classif(out2,newx)
# pred[3,i]=predict.classif(out3,newx)
# pred[4,i]=predict.classif(out4,newx)
# pred[5,i]=predict.classif(out5,newx)
# pred[6,i]=predict.classif(out6,newx)
 pred[7,i]=classif.depth(yy,xx,newx, method = "mode")$group.pred
 pred[8,i]=classif.dist(yy,xx,newx,func=func.trim.mode)$group.pred
 pred[9,i]=classif.dist2(yy,xx,newx,func=func.trim.mode)$group.pred
# if (i>1) print(apply((pred[,1:i]+3)==y[1:i],1,sum))
}

sum((pred[1,]+3)==y)/ng
sum((pred[2,]+3)==y)/ng
sum((pred[3,]+3)==y)/ng
sum((pred[4,]+3)==y)/ng
sum((pred[5,]+3)==y)/ng
sum((pred[6,]+3)==y)/ng
sum((pred[7,]+3)==y)/ng
sum((pred[8,]+3)==y)/ng
sum((pred[9,]+3)==y)/ng
# results
#> sum((pred[1,]+3)==y)/ng
#[1] 0.77
#> sum((pred[2,]+3)==y)/ng
#[1] 0.77
#> sum((pred[3,]+3)==y)/ng
#[1] 0.77
#> sum((pred[4,]+3)==y)/ng
#[1] 0.815
#> sum((pred[5,]+3)==y)/ng
#[1] 0.84
#> sum((pred[6,]+3)==y)/ng
#[1] 0.78
#> sum((pred[7,]+3)==y)/ng
#[1] 0.525
#> sum((pred[8,]+3)==y)/ng
#[1] 0.775
#> sum((pred[9,]+3)==y)/ng
#[1] 0.74
