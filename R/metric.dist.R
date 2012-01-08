metric.dist<-function(x,y=NULL,method="euclidean",p=2,...){
 if (is.vector(x)) x<-matrix(x,nrow=1)

 if (!is.null(y)) {    
    if (is.vector(y)) y<-matrix(y,nrow=1) 
    n<-nrow(y)
    nn<-nrow(x)
    mdist<-as.matrix(dist(rbind(x,y) , method = method, diag = TRUE, upper = TRUE,p=p))[1:nn,(nn+1):(nn+n)] 
    }
 else   mdist<-as.matrix(dist(x, method = method, diag = TRUE, upper = TRUE,p=p))  
   if (is.vector(mdist)) mdist<-matrix(mdist,nrow=1)  
 attr(mdist, "call") <- "metric.dist"
 attr(mdist, "par.metric") <- list(method =method,p=p) 
 return(mdist)
}
