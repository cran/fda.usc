
#library(geoR)
###########################################################################
# Simula una ARMA que permite pasarle model y si el mu=noise y sd=0 
# puede utilizarse para la prediccion
simul.arma<-function(n, model,n.start=1000,mu=0,sd=1) {
    if (!is.list(model)) stop("'model' must be list")
    if (n <= 0L)         stop("'n' must be strictly positive")
    p <- length(model$ar)    
    if (p) {
        minroots <- min(Mod(polyroot(c(1, -model$ar))))
        if (minroots <= 1) 
            stop("'ar' part of model is not stationary")
    }
    q <- length(model$ma)
    if (is.na(n.start)) 
        n.start <- p + q + ifelse(p > 0, ceiling(6/log(minroots)), 
            0)
    if (n.start < p + q) 
        stop("burn-in 'n.start' must be as long as 'ar + ma'")        
    nmax<-n.start+n
    noise <- rnorm(nmax,mu,sd)
    x <- numeric(nmax)  
    maxpq<-max(p,q)
    x[1:(maxpq-1)] <- noise[1:(maxpq-1)]
    for (i in (maxpq+1):nmax)   x[i] <- sum(model$ar * x[(i-1):(i-p)]) + sum(c(1,model$ma)*noise[i:(i-q)])    
    x[(n.start+1):nmax]
    }
# xx<-simul.arma(10,list(ar = c(rho),ma=phi))
# xx<-simul.arma(10,list(ar = c(rho)),1,1,0)
# xx<-simul.arma(10,list(ar = c(rho),ma=.75),2,1,0)

###########################################################################
# Recupera los residuos de un modelo ARIMA
filter.arima<-function(model,x){
##print("entra filter arima")
    p <- length(model$ar)    
    if (p)  x <- filter(x, c(1,-model$ar),method = "convolution",sides = 2L)
    q <- length(model$ma)
    if (q) {
            x[seq_along(model$ma)] <- 0
            x <- filter(x,-model$ma, method = "recursive")
            }
#     x[seq_along(model$ar)] <- 0
##print(x)   
   x #eps[t]<-x[t-1]
}  


###########################################################################
#  rho<-.5
#  nmax<-100   
#  x <- numeric(nmax)
#  noise <-1:nmax
#  x[1] <- noise[1]
#  phi=.75  
#  for (i in 2:nmax)   x[i] <- rho * x[i-1] + noise[i] + phi * noise[i-1]   
#  model<-list(ar = c(rho),ma=phi)
#cbind(x,noise,filter.arima(model,x))
#    if (length(model$ma)) {
#        x <- filter(x, c(1, model$ma), sides = 1L)
#        x[seq_along(model$ma)] <- 0
#    }
# if (length(model$ar)) x <- filter(x, model$ar, method = "recursive")   
###########################################################################


################################################################################
################################################################################
corUnstruc<-function(x){return(x)}
################################################################################
################################################################################
cor.AR<-function(x,order.max=8,p,method="lm"){
# #print("entra cor.AR")
# #print(class(x))
# #print(dim(x))
if (is.vector(x)) {
##print("x is a vector")
 n<-length(x)
 W0<-diag(n)
 x<-x-mean(x)
# #print("AR")
 aa<-ar(x,order.max=order.max, intercept=FALSE,demean=FALSE,method="mle")    
# #print(aa)
 if (aa$order!=0)  {
#  #print("ARIMA")
  aa<-arima(x,order=c(aa$order,0,0),transform.pars=TRUE,include.mean =FALSE,method="CSS")
  
  W0<-toeplitz(ARMAacf(aa$coef,lag.max=c(n-1)))    
#     W0<-toeplitz(ARMAacf(aa$ar,lag.max=n-1))
 }
 else aa$ar<-0
 return(list("W"=W0,"ar"=aa))
 }
# #print("x no is a vector")
# #print(dim(x))
##    W<-solve(W0)  #
#    W <- try(solve(W0),silent=TRUE)
#    if (class(W)=="try-error") {
#      sv<-svd(W0)
#      W<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
#      warning("Inverse of sigma computed by SVD")
#      }
#      #print(round(W0,3)[1:8,1:8])
##   return(W)
if (method=="pc") {
# #print("entra pc")
 x<-t(x)
# x<-cbind(a1,a2,a3) #debe ser un data.frame o una matrix
# meanX<- matrix(apply(x, 2, mean, na.rm = TRUE), ncol = 1)
# Xcen <- sweep(x 2, meanX, FUN = "-")     
 if (!is.matrix(x)) x<-as.matrix(x)
 Xcen<-x # si son los residuos ya vienen centrados (pero no por grupos)
 eigenres <- svd(Xcen)
    v <- eigenres$v
    u <- eigenres$u
    d <- eigenres$d
    l<-1
    lenl<-length(l)
    pc.fdata <- u[, l, drop = FALSE] %*% (diag(lenl) * d[l]) %*% v[l, , drop = FALSE]
    aa<-ar(pc.fdata[,1],order.max=order.max)
##print(aa)
    n<-length(pc.fdata[,1])
    W0<-diag(n)
    if (aa$order!=0)  W0<-toeplitz(ARMAacf(aa$ar,lag.max=n-1))
    else aa$ar<-0
# #print("sale pc")
# #print(dim(W0))
   return(list("W"=W0,"ar"=aa,"pc"= eigenres ))
 }
if (method=="lm") {
# #print("entra lm")
 x<-t(x)
# estimar un lm  de orden?
 n<-nrow(x)
 np<-ncol(x)
# #print(missing(p))
if (missing(p)) p<-1:order.max
else {
      order.max<-p
      p<-1:p      
      }      
 lenp<-order.max
 ##print(p)
 maxp<-max(p)
 aabb<-list()
 nam<-paste("x.",p,sep="")
# aic<-Inf
aabb<-list()
 aabb[["y"]]<-as.vector(x[(order.max+1):n,])
 best<-lm(y~1,data=aabb)
 aic<-AIC(best)
 bestp<-0
 coeff<-0
#   #print(aic)
 for (i in 1:lenp) {
 ##print("entra iiiiii")
   aabb[[nam[i]]]<-as.vector(x[(order.max+1-i):(n-i),])
#   #print("oooo")
   f<-as.formula(paste("y~-1+",paste(nam[1:i],collapse="+"),sep=""))
#  #print(f)
   aa<-lm(f,data=aabb)
   aic2<-AIC(aa)
   if (aic2<aic) {
                  bestp<-i
                  best<-aa
                  aic<-aic2
                  }
 }
 if (bestp!=0) {
 coeff <- best$coefficients
# #print(coeff)
 best$res.lm<-aabb[1:(bestp+1)]
 best$res.x<-x[n:(n-bestp+1),,drop=FALSE]
 W0<-toeplitz(ARMAacf(coeff,lag.max=n-1))
 }
 else  {
    W0<-diag(np)
    coeff<-0
    best$res.x<-matrix(0,ncol=ncol(x))
    }

##print(best)
# #print(dim(W0))
#    #print("sale lm")
 return(list("W"=W0,"lm"=best))
 }
}

################################################################################
################################################################################
cor.ARMA<-function(x,p,d=0,q=0,method="lm",order.max=1){
# #print("entra cor.ARMA")
# #print(class(x))
# #print(dim(x))

if (is.vector(x)) {
##print("x is a vector")
 n<-length(x)
 W0<-diag(n)
 if (missing(p)) p<-1
 aa<-arima(x,order=c(p,d,q),include.mean =FALSE,transform.pars=TRUE)
     W0<-toeplitz(ARMAacf(aa$coef,lag.max=c(n-1)))
# #print("sale AR")
 return(list("W"=W0,"ar"=aa))
}
# #print("x no is a vector")
# #print(dim(x))
##    W<-solve(W0)  #
#    W <- try(solve(W0),silent=TRUE)
#    if (class(W)=="try-error") {
#      sv<-svd(W0)
#      W<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
#      warning("Inverse of sigma computed by SVD")
#      }
#      #print(round(W0,3)[1:8,1:8])
##   return(W)
if (method=="pc") {
# #print("entra pc")
 x<-t(x)
# x<-cbind(a1,a2,a3) #debe ser un data.frame o una matrix
# meanX<- matrix(apply(x, 2, mean, na.rm = TRUE), ncol = 1)
# Xcen <- sweep(x 2, meanX, FUN = "-")     
 if (!is.matrix(x)) x<-as.matrix(x)
 Xcen<-x # si son los residuos ya vienen centrados (pero no por grupos)
 eigenres <- svd(Xcen)
    v <- eigenres$v
    u <- eigenres$u
    d <- eigenres$d
    l<-1
    lenl<-length(l)
    pc.fdata <- u[, l, drop = FALSE] %*% (diag(lenl) * d[l]) %*% v[l, , drop = FALSE]
    aa<-ar(pc.fdata[,1],order.max=order.max)
##print(aa)
    n<-length(pc.fdata[,1])
    W0<-diag(n)
    if (aa$order!=0)  W0<-toeplitz(ARMAacf(aa$ar,lag.max=n-1))
    else aa$ar<-0
   return(list("W"=W0,"ar"=aa,"pc"= eigenres ))
 }
if (method=="lm") {
# #print("entra lm")
 x<-(x)
# estimar un lm  de orden?
 n<-nrow(x)
 np<-ncol(x)
# #print(missing(p))
if (missing(p)) p<-1:order.max
else {
      order.max<-p
      p<-1:p      
      }      
 lenp<-order.max
 maxp<-max(p)
 aabb<-list()

 nam<-paste("x.",p,sep="")
 aabb[["y"]]<-as.vector(x[(maxp+1):n,])
 for (i in 1:maxp) aabb[[nam[i]]]<-as.vector(x[(maxp+1-i):(n-i),])
#   #print("oooo")
   f<-as.formula(paste("y~-1+",paste(nam[p],collapse="+"),sep=""))
   aa<-lm(f,data=aabb)
   coeff <- aa$coefficients
 aa$res.lm<-aabb[1:(maxp+1)]
 aa$res.x<-x[n:(n-maxp+1),,drop=FALSE]  #la matriz de los ultimos residuo
 W0<-toeplitz(ARMAacf(coeff,lag.max=n-1))
# #print("sale lm")
 return(list("W"=W0,"lm"=aa))
 }
}
#######################################################################


###############################################################################
###############################################################################
corExpo<-function(xy,range, method = "euclidean",p=2){
 if (is.data.frame(xy)) xy<-as.matrix(xy)
 if (class(xy)=="dist") dxy<-as.matrix(xy)
 else dxy<- as.matrix(dist(xy,method=method,p=p,diag =TRUE, upper = TRUE))
# como afecta los residuos?
vdxy<-as.vector(dxy)
if (missing(range)) range<-quantile(vdxy,.9)/3
else range/3
  W0=exp(-dxy/range)
  return(W0)
  n<-nrow(xy)
##   W<-solve(W0)
#    W <- try(solve(W0),silent=TRUE)
#    if (class(W)=="try-error") {
#      sv<-svd(W0)
#      W<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
#      warning("Inverse of sigma computed by SVD")
#      }
#  #print(round(W0,3)[1:3,1:3])
##   return(W)   
 }


# 2016/10/31 se elimina cor.Exp y corCloud 
###############################################################################
