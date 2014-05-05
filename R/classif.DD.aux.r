################################################################################
# Auxiliary functions for classif.DD
# File created by Manuel Oviedo de la Fuente  using code from paper:
# Li, J., P.C., Cuesta-Albertos, J.A. and Liu, R. 
# DD--Classifier: Nonparametric Classification Procedure Based on DD-plot. 
# Journal of the American Statistical Association (2012), Vol. 107, 737--753. 
##########################################################################################################################################
RR <- function(x,a){
   y <- 0
   kk <- length(a)
   for (i in 1:kk){
        y <- y+ a[i]*(x^i)
   }
   return(y)
}

AMCR0 <- function(a,dep,ind,tt){
  p<-sum(ind[,1])/nrow(ind)
  amcr <- p*mean(1/(1+exp(-tt*(dep[ind[,1],2]-sapply(dep[ind[,1],1],RR,a=a)))))+
        (1-p)*mean(1/(1+exp(-tt*(sapply(dep[ind[,2],1],RR,a=a)-dep[ind[,2],2]))))
  return(amcr)
}


quad.fit0 <- function(index,dep,ind){
   Df<-dep[,1]
   Dg<-dep[,2]
    A <- matrix(c(Df[index[1]],(Df[index[1]])^2, Df[index[2]],(Df[index[2]])^2),byrow=TRUE,2,2)
     w <- c(Dg[index[1]],Dg[index[2]])
    a <- try(solve(A,w),silent=TRUE)
    if (class(a)=="try-error") {
     sv<-svd(A)
     a<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u))%*%w)
#    warning("Inverse of sigma computed by SVD")
    }
    mcr<-MCR0.p(a,dep,ind)
    return(mcr)
}

cubic.fit0 <- function(index,dep,ind){
   Df<-dep[,1]
   Dg<-dep[,2]
    A <- matrix(c(Df[index[1]],(Df[index[1]])^2,(Df[index[1]])^3,
          Df[index[2]],(Df[index[2]])^2,(Df[index[2]])^3,
          Df[index[3]],(Df[index[3]])^2,(Df[index[3]])^3),byrow=TRUE,3,3)
    w <- c(Dg[index[1]],Dg[index[2]],Dg[index[3]])
    a <- try(solve(A,w),silent=TRUE)

    if (class(a)=="try-error") {
     sv<-svd(A)
     a<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u))%*%w)
 #    warning("Inverse of sigma computed by SVD")
    }
    mcr<-MCR0.p(a,dep,ind)
    return(mcr)
}

##########################################################################################################################################
# the function calculates the empirical misclassification rate from xx and yy based on the polynomial curve with coefficients being a in the DD-plot
# (the points on the curve are classified using KNN) NO

MCR0.p <- function(a,dep,ind){
 dm<-dim(ind)
 mis<-0
 p<-colMeans(ind)
 ahorro<-sapply(dep[ind[,1],1],RR,a=a)
 ahorro2<-sapply(dep[ind[,2],1],RR,a=a)
 for (i in 1:(dm[2]-1)) {
  mis<-p[i]*mean(ahorro<dep[ind[,1],2])+(1-p[i])*mean(ahorro2>dep[ind[,2],2])
  }
 mis
 }


##########################################################################################################################################
# the function calculates the empirical misclassification rate from xx and yy based on the linear line with slope being k in the DD-plot
# (the points on the linear line are classified using KNN) NO

MCR0 <- function(k,dep,ind){
dm<-dim(ind)
mis<-0
p<-colMeans(ind)
for (i in 1:(dm[2]-1)) {
  mis<-mis+p[i]*mean((k*dep[ind[,i],1])<dep[ind[,i],2])+(1-p[i])*mean((k*dep[-ind[,i],1])>dep[-ind[,i],2])
  }
 mis
 }
##########################################################################################################################################
##########################################################################################################################################