################################################################################
# File created by Manuel Oviedo de la Fuente  using code from paper:
# Li, J., P.C., Cuesta-Albertos, J.A. and Liu, R. 
# DD--Classifier: Nonparametric Classification Procedure Based on DD-plot. 
# Journal of the American Statistical Association (2012), Vol. 107, 737--753. 
################################################################################
        
###############################################################################################################
#mdepth.HS:  calculates the half-space depth (HS) of the points in x w.r.t. xx based on projections 
        #xx is a d-dimension multivariate sample, a d-column matrix
        #x is a set of points, a d-column matrix, x can be missing
        #proj are the directions for projections    
mdepth.HS <-function(x, xx=x,proj=50,scale=FALSE,xeps=1e-15,random=FALSE)
{
      if (is.vector(x)) {
         x<-matrix(x,nrow=1)
#         ans<-pmin(sum(x<=xx)/m,sum(x>=xx)/m)
         Fn=ecdf(xx)
         ans=pmin(Fn(x),(1-Fn(x-xeps)))
         if (scale) ans<-ans*2   
         return(invisible(list(dep = ans, Fn=Fn)))   
        }
       if ( is.null(rownames(x)))  rownames(x)<-1:nrow(x)
       nms<-rownames(x)
       m0<-nrow(x)
       xx<-na.omit(xx)
       x<-na.omit(x)
       nas<-na.action(x)
       nullans<-!is.null(nas)        
        n <- nrow(x)
        d <- dim(xx)[2]
        lenn<-length(proj)
        if (lenn==1) {
          mm<-proj[1]
          if (d==2 & !random) {
              sq<-seq(0,2*pi,len=mm)
#              proj<-as.matrix(expand.grid(cos(sq),sin(sq)))             
               proj2d<-function(angl){matrix(c(cos(angl),sin(angl)),2)}
               proj<-t(sapply(sq,proj2d   ))
             }
          else {
          if (d==3 & !random) {
			      mmr=floor(sqrt(mm))+1
      		  phi=seq(0,2*pi,len=mmr)
			      theta=seq(0,pi,len=mmr)
            exgrid=expand.grid(phi=phi,theta=theta)
			      proj=cbind(sin(exgrid$theta)*cos(exgrid$phi),sin(exgrid$theta)*sin(exgrid$phi),cos(exgrid$theta))
            mm<-nrow(proj)          
            }
          else{             
          warning("Method based on Random Projections")
          u <- matrix(runif(d*mm,-1,1),mm,d)
          norm <- sqrt(rowSums(u*u))
          proj <- (u/norm) 
         }
        }            
        }
        else   mm<-nrow(proj)
      out <- matrix(0, mm,n) 
      if (d==3 & !random) {
               for(i in 1:mm) {
#                 z<-t(a[,,i]%*%t(x))# este calculo esta malm  debe dar nproj*ndatos
#                 z2<-t(a[,,i]%*%t(xx))
        z  <- proj %*% t(x)            
        z2 <- proj %*% t(xx)
                  
                 out[i,] <- sapply(z[i,], function(y) min(sum(y<=z2[i,])/n,sum(y>=z2[i,])/n))
               }       
      }
      else{     
        z  <- proj %*% t(x)            
        z2 <- proj %*% t(xx)
        for(i in 1:mm) {
          out[i,] <- sapply(z[i,], function(y) min(sum(y<=z2[i,])/n,sum(y>=z2[i,])/n))
        }        
        }
    ans = as.vector(apply(out,2,min))        
   if (scale)          ans<-ans*2
  if  (nullans){
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-ans 
        ans<-ans1      
        }
   names(ans)<-nms   
   out <- list(dep = ans, proj = proj, x=x,xx=xx,name="HS")
   return(invisible(out))
}

#names
# revisar 3d sapply
#  a<-mdepth.HS(iris[,1:2],random=F)
# b<-mdepth.HS(iris[,1:2],random=T)  
#  plot(a$dep,b$dep)
# cor(a$dep,b$dep)

#   a<-mdepth.HS(iris[,1:3],random=F,proj=5)
#   b<-mdepth.HS(iris[,1:3],random=T)  
#   plot(a$dep,b$dep)
#  cor(a$dep,b$dep)


#################################################################################################################
# depth.RD: calculates the projection depth (PD) of the points in x w.r.t. xx based on random projections proj
        #xx is a d-dimension multivariate sample, a d-column matrix
        #x is a set of points, a d-column matrix, x can be missing
        #proj are the directions fo random projections
        #trim the alpha of the trimming
        #draw=TRUE, draw the points in a gray scale of its depth, the sample median (in red) and trimmed mean (in yellow)
     
mdepth.RP<-function(x, xx=x,proj=50,scale=FALSE){

 if (is.vector(x)){
   if (all(xx==x)) x<-xx<-matrix(x,ncol=1)#stop("One of x or xx must be a matrix")
   else   {
 	m2=ncol(xx)
	if (length(x)!=m2) stop("Length of x does not match with dimension of xx")
	x=matrix(x,ncol=m2)
   }
 }
     if ( is.null(rownames(x)))  rownames(x)<-1:nrow(x)
       nms<-rownames(x)
       m0<-nrow(x)
       xx<-na.omit(xx)
       x<-na.omit(x)
       nas<-na.action(x)
       nullans<-!is.null(nas) 
        n <- nrow(x)
        d <- ncol(x)
        lenn<-length(proj)
        if (lenn==1) {
          u <- matrix(runif(d*proj,-1,1),proj,d)
          norm <- sqrt(rowSums(u*u))
          proj <- u/norm
        }
        z <- proj %*% t(xx)
        z1 <- proj %*% t(x)        
        mm<-nrow(proj)
		pdep=matrix(NA,nrow=n,ncol=mm)
		for (i in 1:mm){
		pdep[,i]=mdepth.TD1(z1[i,],z[i,],scale=scale)$dep
#        m1 <- m2 <- rep(0, mm)
#        for(i in 1:mm) {
#				 m1[i]= min(z[i,])
#				 m2[i]= max(z[i,])
#                m1[i] <- median(z[i,  ])
#                m2[i] <- median(abs(z[i,  ] - m1[i]))
        }
        out1 <- apply(pdep,1,mean)  
#        for(j in 1:n) {  out1[j] <- max(abs(z1[, j] - m1)/m2)       }             
#        for(j in 1:n) {  out1[j] <- mean(abs(z1[, j] - m1)/(m2-m1))       }             
		ans=out1
#        ans = 1/(1 + out1)        
#        if (scale){ 
#          ans<-ans/max(ans) #*2 =/.5
#        }   
  if  (nullans){
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-ans 
        ans<-ans1
        }
     names(ans)<-nms    
     out <- list( dep = ans,  proj = proj,x=x,xx=xx,name="RP")
    return(invisible(out))
}

#################################################################################
#################################################################################
#mdepth.MhD:  calculates the Mahalanobis depth (MhD) of the points in x w.r.t. xx
        #xx is a d-dimension multivariate sample, a d-column matrix
        #x is a set of points, a d-column matrix, x can be missing
        #trim the alpha of the trimming
        #draw=TRUE, draw the points in a gray scale of its depth, the sample median (in red) and trimmed mean (in yellow)
       
mdepth.MhD <- function(x,xx=x,scale=FALSE){
 m0 <- nrow(x)
 if (!is.vector(x) & is.null(rownames(x)))  rownames(x)<-1:nrow(x)
 nms<-rownames(x)
 x<-na.omit(x)    
 xx<-na.omit(xx)
 nas<-na.action(x)
 nullans<-!is.null(nas) 
 if (is.vector(x)) {
              D<- (1+(x-(mean(xx)))^2/sd(xx)^2)
              }
 else{

  n <- nrow(xx)
	m <- nrow(x)
	d<-ncol(x)
	mu <-colMeans(xx)
	sigma <- cov(xx)
	D <-rep(0,m)
  sigma.inv <- try(solve(sigma),silent=TRUE)#new
  
  if (!is.matrix(sigma.inv)) {
     sv<-svd(sigma)    
     sigma.inv<-sv$v%*%diag(1/sv$d)%*%t(sv$u)
     warning("Inverse of sigma computed by SVD")
    }
   D <- 1+apply(t(x)-mu,2, function(x) t(x)%*%sigma.inv%*%x)
   }
ans<-1/D      
if  (nullans){
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-1/D 
        ans<-ans1
        }
  names(ans)<-nms   
#   if (scale) {        ans <- ans/max(ans)    }
  out <- list( dep = ans,  x=x,xx=xx,name="MhD")
   return(invisible(out))
}  

  
#################################################################################
#################################################################################
