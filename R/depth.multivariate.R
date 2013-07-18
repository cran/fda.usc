################################################################################
# File created by Manuel Oviedo de la Fuente  using code from paper:
# Li, J., P.C., Cuesta-Albertos, J.A. and Liu, R. 
# DD--Classifier: Nonparametric Classification Procedure Based on DD-plot. 
# Journal of the American Statistical Association (2012), Vol. 107, 737--753. 
################################################################################
        
#################################################################################################################
#depth.MhD:  calculates the Mahalanobis depth (MhD) of the points in x w.r.t. xx
        #xx is a d-dimension multivariate sample, a d-column matrix
        #x is a set of points, a d-column matrix, x can be missing
        #trim the alpha of the trimming
        #draw=TRUE, draw the points in a gray scale of its depth, the sample median (in red) and trimmed mean (in yellow)
       
depth.MhD <- function(x,xx=x,trim=0.25,draw=FALSE){
 if (is.vector(x)) x<-matrix(x,nrow=1)
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
   D2 <- 1+apply(t(xx)-mu,2, function(x) t(x)%*%sigma.inv%*%x)
   ans <- 1/D
   ans2 <- 1/D2
   k=which.max(ans2)
   med=xx[k,]
   lista=which(ans2>=quantile(ans2,probs=trim,na.rm=TRUE))
       if (length(lista)==1) {
            mtrim<-xx[lista,]
            if (draw) {draw=FALSE;warning("The plot is not shown")}
        }
        else mtrim=apply(xx[lista,],2,mean,na.rm=TRUE)
      dd=0      
        if (draw==TRUE) {
          if (d==2) dd=2
          if (d>3)  dd=3
          if (d>5)  dd=4          
          }        
        else {
          switch(draw,"pairs"={dd=3},"stars"={dd=4})       
        }   
        if (dd!=0) {
         tr <- paste("PD.trim", trim * 100, "%", sep = "")
         ind1 <- !is.na(ans2)
         cgray = 1 - (ans2 - min(ans2, na.rm = TRUE))/(max(ans2,
             na.rm = TRUE) - min(ans2, na.rm = TRUE))
         if (is.data.frame(x)) nam<-names(xx)
         if (is.matrix(xx)) nam<-colnames(xx)
         if (dd==2) {
          plot(xx[ind1, 1],xx[ind1, 2], col = gray(cgray[ind1]),main = "Projection Depth",xlab=nam[1],ylab=nam[2])
          points(mtrim[1],mtrim[2], lwd = 2, col = "yellow",pch=16)
          points(med[1],med[2], col = "red", lwd = 2,pch=17)
          legend("topleft", legend = c(tr, "Median"), pch=16:17, box.col=0,col = c("yellow", "red")) }
          if (dd==3) pairs(xx[ind1, ],pch=1,col = gray(cgray[ind1]), main = "Projection Depth") 
          if (dd==4) stars(xx[ind1, ],col.stars = gray(cgray[ind1]))
    }
    return(invisible(list(median = med, lmed = k, mtrim =matrix(mtrim,nrow=1),ltrim = lista, dep = ans, dep.ori = ans2)))
}   





#################################################################################################################
#depth.HD:  calculates the half-space depth (HD) of the points in x w.r.t. xx based on random projections proj
        #xx is a d-dimension multivariate sample, a d-column matrix
        #x is a set of points, a d-column matrix, x can be missing
        #proj are the directions fo random projections
        #trim the alpha of the trimming
        #draw=TRUE, draw the points in a gray scale of its depth, the sample median (in red) and trimmed mean (in yellow)
      
depth.HD <-function(x, xx=x,trim=0.25,proj=NULL,draw=FALSE)
{
        if (is.vector(x)) x<-matrix(x,nrow=1)
        n <- nrow(x)
        d <- dim(xx)[2]
        m <- length(xx[, 1])
        if (is.null(proj)) {
          u <- matrix(runif(d*500,-1,1),500,d)
          norm <- sqrt(rowSums(u*u))
          proj <- u/norm
        }
        mm<-nrow(proj)
        z <- proj %*% t(xx)
        z2 <- proj %*% t(x)
        out <- matrix(0, mm,m)
        out2 <- matrix(0, mm,n)
        for(i in 1:mm) {
        out[i,] <- sapply(z[i,], function(y) min(sum(y<=z[i,])/m,sum(y>=z[i,])/m))
        out2[i,] <- sapply(z2[i,], function(y) min(sum(y<=z[i,])/m,sum(y>=z[i,])/m))
        }
        ans = as.vector(apply(out,2,min))
        ans2 = as.vector(apply(out2,2,min))        
        k=which.max(ans)
        med=xx[k,]
        lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
        if (length(lista)==1) {
            mtrim<-xx[lista,]
            if (draw) {draw=FALSE;warning("The plot is not shown")}
        }
        else mtrim=apply(xx[lista,],2,mean,na.rm=TRUE)
        dd=0      
        if (draw==TRUE) {
          if (d==2) dd=2
          if (d>3)  dd=3
          if (d>5)  dd=4          
          }        
        else {
          switch(draw,"pairs"={dd=3},"stars"={dd=4})       
        }   
        if (dd!=0) {
         tr <- paste("PD.trim", trim * 100, "%", sep = "")
         ind1 <- !is.na(ans)
         cgray = 1 - (ans - min(ans, na.rm = TRUE))/(max(ans,
             na.rm = TRUE) - min(ans, na.rm = TRUE))
         if (is.data.frame(x)) nam<-names(xx)
         if (is.matrix(xx)) nam<-colnames(xx)
         if (dd==2) {
          plot(x[ind1, 1],xx[ind1, 2], col = gray(cgray[ind1]),main = "Projection Depth",xlab=nam[1],ylab=nam[2])
          points(mtrim[1],mtrim[2], lwd = 2, col = "yellow",pch=16)
          points(med[1],med[2], col = "red", lwd = 2,pch=17)
          legend("topleft", legend = c(tr, "Median"), pch=16:17, box.col=0,col = c("yellow", "red")) }
          if (dd==3) pairs(xx[ind1, ],pch=1,col = gray(cgray[ind1]), main = "Projection Depth") 
          if (dd==4) stars(xx[ind1, ],col.stars = gray(cgray[ind1]))
    }
    return(invisible(list(median = med, lmed = k, mtrim = matrix(mtrim,nrow=1),
        ltrim = lista, dep = ans2,dep.ori = ans, proj = proj)))
}


#################################################################################################################
# depth.PD: calculates the projection depth (PD) of the points in x w.r.t. xx based on random projections proj
        #xx is a d-dimension multivariate sample, a d-column matrix
        #x is a set of points, a d-column matrix, x can be missing
        #proj are the directions fo random projections
        #trim the alpha of the trimming
        #draw=TRUE, draw the points in a gray scale of its depth, the sample median (in red) and trimmed mean (in yellow)
     
depth.PD <-function(x, xx=x,trim=0.25,proj=NULL,draw=FALSE){
        if (is.vector(x)) x<-matrix(x,nrow=1)
        n <- nrow(x)
        d <- ncol(x)
        if (is.null(proj)) {
          u <- matrix(runif(d*500,-1,1),500,d)
          norm <- sqrt(rowSums(u*u))
          proj <- u/norm
        }
        z <- proj %*% t(xx)
        mm<-nrow(proj)
        m1 <- m2 <- rep(0, mm)
        for(i in 1:mm) {
                m1[i] <- median(z[i,  ])
                m2[i] <- median(abs(z[i,  ] - m1[i]))
        }
        m <- length(xx[, 1])
        z <- proj %*% t(xx)
        z1 <- proj %*% t(x)        
        out <- rep(0, m)
        out1 <- rep(0, n)  
        for(j in 1:m) {  out[j] <- max(abs(z[, j] - m1)/m2)       }
        for(j in 1:n) {  out1[j] <- max(abs(z1[, j] - m1)/m2)       }             
        ans = 1/(1 + out)
        ans1 = 1/(1 + out1)        
        k=which.max(ans)
        med=xx[k,]
        lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
        if (length(lista)==1) {
          mtrim<-xx[lista,]
          if (draw) {draw=FALSE;warning("The plot is not shown")}
        }
        else mtrim=apply(xx[lista,],2,mean,na.rm=TRUE)
        dd=0      
        if (draw==TRUE) {
          if (d==2) dd=2
          if (d>3)  dd=3
          if (d>5)  dd=4          
          }        
        else {
          switch(draw,"pairs"={dd=3},"stars"={dd=4})       
        }   
        if (dd!=0) {
         tr <- paste("PD.trim", trim * 100, "%", sep = "")
         ind1 <- !is.na(ans)
         cgray = 1 - (ans - min(ans, na.rm = TRUE))/(max(ans,
             na.rm = TRUE) - min(ans, na.rm = TRUE))
         if (is.data.frame(xx)) nam<-names(xx)
         if (is.matrix(xx)) nam<-colnames(xx)
         if (dd==2) {
          plot(xx[ind1, 1],xx[ind1, 2], col = gray(cgray[ind1]),main = "Projection Depth",xlab=nam[1],ylab=nam[2])
          points(mtrim[1],mtrim[2], lwd = 2, col = "yellow",pch=16)
          points(med[1],med[2], col = "red", lwd = 2,pch=17)
          legend("topleft", legend = c(tr, "Median"), pch=16:17, box.col=0,col = c("yellow", "red")) }
          if (dd==3) pairs(x[ind1, ],pch=1,col = gray(cgray[ind1]), main = "Projection Depth") 
          if (dd==4) stars(x[ind1, ],col.stars = gray(cgray[ind1]))
    }
    return(invisible(list(median = med, lmed = k, mtrim = matrix(mtrim,nrow=1),
        ltrim = lista, dep = ans1, dep.ori = ans, proj = proj)))
}

