metric.dist<-function(x,y=NULL,method="euclidean",p=2,...){
 if (!is.null(y)) {    
    n<-nrow(y)
    nn<-nrow(x)
    mdist<-as.matrix(dist(rbind(x,y) , method = method, diag = TRUE, upper = TRUE,p=p))[1:nn,(nn+1):(nn+n)] 
    }
 else   mdist<-as.matrix(dist(x, method = method, diag = TRUE, upper = TRUE,p=p))  
 attr(mdist, "call") <- "metric.dist"
 attr(mdist, "par.metric") <- list(method =method,p=p) 
 return(mdist)
}