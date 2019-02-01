depth.RT  <-function (fdataobj,fdataori=fdataobj, trim = 0.25, nproj = 10, proj = 1, xeps = 1e-07, 
    draw = FALSE, ...) 
{
    if (!is.fdata(fdataobj)) fdataobj = fdata(fdataobj)
    if (!is.fdata(fdataori)) fdataobj=fdata(fdataori)    
	if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
 nms<-rownames(fdataobj$data)
 m0<-nrow(fdataobj)
 fdataobj2<-fdataobj
 fdataobj<-na.omit.fdata(fdataobj)
 nas<-na.action(fdataobj)
 nullans<-!is.null(nas) 
 
 fdataori=na.omit.fdata(fdataori)
 data  <- fdataobj[["data"]]
    data2 <- fdataori[["data"]]    
    n <- nrow(fdataobj)
    m <- ncol(fdataobj)
    m2<-ncol(fdataori)
    n2<-nrow(fdataori)
    if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS OF fdataobj")
    if (is.null(n2) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS fdataori")
    tt = fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names1 <- names2 <- names <- fdataobj[["names"]]
    names1$main <- "depth.RT median"
    names2$main <- paste("RT trim", trim * 100, "%", sep = "")
   
    if (is.fdata(proj)) {
        nproj <- nrow(proj)
        if (fdataobj$argvals != proj$argvals || ncol(fdataobj) != 
            ncol(proj)) 
            stop("Error en proj dimension")
        z <- proj
    }
    else {
        z <- rproc2fdata(nproj, tt, sigma = proj, norm = TRUE,...)
    }
    Fn <- list()
	Prod=inprod.fdata(fdataobj,z)
	Prod2=inprod.fdata(fdataori,z)
    dep = array(NA, dim = c(nproj, n))
  
    for (j in 1:nproj) {
        Fn[[j]] = ecdf(Prod2[,j])  
        dep[j, ] = pmin(Fn[[j]](Prod[,j]-xeps) ,(1 - Fn[[j]](Prod[,j]-xeps)))       
    }
#    print(dep)
    dep = apply(dep, 2, min)
    if (nullans) {
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-dep
        dep<-ans1      
        }
    names(dep)<-nms      
    k = which.max(dep)
    med = fdataobj2[k]
    nl = length(trim)
	lista=vector("list",nl)
    tr <- paste("RT.tr", round(trim * 100,2), "%", sep = "")	
	if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
    mtrim = fdata(matrix(NA, nrow = nl, ncol = m),tt,rtt,names)
    for (j in 1:length(trim)) {
        lista[[j]] = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
        if (length(lista[[j]])==1) {
          mtrim$data[j,]<-fdataobj2[lista[[j]]]$data
          if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
        }
        else mtrim$data[j,]=func.mean(fdataobj2[lista[[j]]])$data       
    }
    rownames(med$data) <- "RT.med"
    rownames(mtrim$data) <- tr
	out=list(median = med, lmed = k, mtrim = mtrim, 
        ltrim = if (nl==1) unlist(lista) else lista, dep = dep, proj = z, Fn = Fn,name="RT")
	out$trim <- trim
	out$name <- "RT"
	out$fdataobj <- fdataobj
	out$fdataori <- fdataori
	class(out)="depth"
    if (draw) {
		plot.depth(out)
    }
 return(invisible(out))
}
