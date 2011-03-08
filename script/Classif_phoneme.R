################################################################################
# We present three non parametric methods for discriminating functional data
# 1.k-Nearest Neighbor method for fucntional data
# 2.Non parametric kernel method for fucntional data
# 3.Non parametric kernel method for basis representation of fucntional data
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
out1=classif.knn.fd(mlearn,glearn,knn=k)

#2. the argument "h" is vector of bandwidth considered
h=h.default(mlearn,prob=c(0.001,0.1),len=10)
out2=classif.kernel.fd(mlearn,glearn,h=h)

#3. the argument "type.basis" is the name of the basis
# where the functional data is represented
out3=classif.kernel.fb(mlearn,glearn,type.basis="fourier",h=h)

################################################################################
# Summarizes information from classification methods of class "classif.fd", shows:
# -Probability of correct classification by group: prob.classification.
# -Confusion matrix between the theoretical groups and estimated groups.
# -Highest probability of correct classification max.prob.
################################################################################
#1. If the object is returned from the function: classif.knn.fd()
# -Vector of probability of correct classification by number of neighbors knn.
# -Optimal number of neighbors: knn.opt.
summary.classif.fd(out1)

#2. If the object is returned from the function: classif.kernel.fd
#   -Vector of probability of correct classification by banwidth h.
#  -Functional measure of closeness (optimal distance, h.opt).
summary.classif.fd(out2)

#3. If the object is returned from the function: classif.kernel.fb
#   -Matrix of probability of correct classification by number of basis num.basis and bandwidth h.
#   -Functional measure of proximity(optimal distance, h.opt).
#   -Optimal number of basis: basis.opt and optimal bandwidth: h.opt.
summary.classif.fd(out3)


################################################################################
# Read test curves
mtest<-phoneme[["test"]]

#	PREDICTION  classmember from new curves
pred1=predict.classif.fd(out1,mtest,TRUE)
names(pred1)
#estimated groups: pred1$group.pred
#for each curve shows the probability of each group: pred1$prob.group
pred2=predict.classif.fd(out2,mtest,TRUE)
pred3=predict.classif.fd(out3,mtest,TRUE)

# Read theoretical groups
gtest<-phoneme[["classtest"]]

# Confusion matrix between the theoretical groups and estimated groups
t1<-table(pred1$group.pred,gtest)
t2<-table(pred2$group.pred,gtest)
t3<-table(pred3$group.pred,gtest)

# Missclassification  rate
ntest<-length(gtest)
MisclasPredict1 <- sum(pred1$group.pred != gtest)/ntest
MisclasPredict2 <- sum(pred2$group.pred != gtest)/ntest
MisclasPredict3 <- sum(pred3$group.pred != gtest)/ntest
c(MisclasPredict1,MisclasPredict2,MisclasPredict3)

# Probability of correct classification by group: prob.classification.
prob.class<-cbind(diag(t1),diag(t2),diag(t2))/50
#rownames(prob.class)<-c("aa","ao","dcl","iy","sh")
colnames(prob.class)<-c("knn.fd","kernel.fd","kernel.fb")
prob.class
matplot(prob.class,ylab="Probability Correct Classification",xlab="Phonemes")



