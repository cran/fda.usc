#######################################################################
rproc2fdata=function(n,t,mu=rep(0,length(t)),sigma=1,
par.list=list("scale"=1,"theta"=1/3),norm=FALSE,...) {
p=length(t)
sigma2<-sigma
if (is.null(par.list$mu0)) par.list$mu0<-0
if (is.null(par.list$theta)) par.list$theta<-1/3
if (is.null(par.list$scale)) par.list$scale<-1
if (is.character(sigma)) {
 type.proc<-c("brownian","wiener","OU","OrnsteinUhlenbeck","vexponential")
 if (!is.element(sigma,type.proc)) stop("Error in sigma argument label")
 ss<-which(sigma==type.proc)
  sigma=switch(sigma,"brownian"=par.list$scale*outer(t,t,function(u,v){pmin(u,v)}),
         "wiener"= par.list$scale*outer(t,t,function(u,v){pmin(u,v)}),
         "OU"= {m<-100;t2<-t+m
         par.list$scale/(2*par.list$theta)*outer(t2,t2,function(u,v){
            exp(-par.list$theta*(u+v))*(exp(2*par.list$theta*pmin(u,v))-1)})},
"OrnsteinUhlenbeck"={m<-100;t2<-t+m
            par.list$scale/(2*par.list$theta)*outer(t2,t2,function(u,v){
            exp(-par.list$theta*(u+v))*(exp(2*par.list$theta*pmin(u,v))-1)})},
             "vexponential"=par.list$theta*outer(t,t,function(u,v){(1-exp(-3*abs(u-v)/par.list$scale))})
             )
  sigma<-t(sigma)

}
#   else if (sigma=="vgaussian") {
#   sigma=par.list$theta*outer(t,t,function(u,v){(1-exp(-3*(u-v)^2/par.list$scale))})}
else {
 if   (is.matrix(sigma)) {
  if (dim(sigma)[2]!=p) stop("Error in sigma argument")
 }
 else if (length(sigma)==1) sigma<-diag(p)*sigma
 else stop("Error in sigma argument")
}
C=svd(t(sigma))
L=C$u%*%diag(sqrt(C$d))
X=matrix(rnorm(n*p),ncol=p)
X=t(L%*%t(X))
X=sweep(X,2,mu,"+")
X=fdata(X,t)
if (norm) {
if (sigma2[1]=="brownian") print("The normalization is not done")
else{
        no <- norm.fdata(X,...)
        X$data <- sweep(X$data, 1, drop(no), "/")
        }
}
return(X)
}

