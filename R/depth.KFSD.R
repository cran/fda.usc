depth.KFSD=function (fdataobj, fdataori = fdataobj, trim = 0.25,
                     h=NULL,scale = FALSE, draw = FALSE){
    if (is.fdata(fdataobj) & is.fdata(fdataori)) {
        if (is.null(rownames(fdataobj$data))) 
            rownames(fdataobj$data) <- 1:nrow(fdataobj$data)
        nms <- rownames(fdataobj$data)
        m0 <- nrow(fdataobj)
        fdataobj <- na.omit.fdata(fdataobj)
        fdataori <- na.omit.fdata(fdataori)
        nas <- na.action(fdataobj)
        nullans <- !is.null(nas)
        names1 <- names2 <- names <- fdataobj[["names"]]
        names1$main <- "depth.KFSD median"
        names2$main <- paste("depth.KFSD trim ", trim * 100, 
            "%", sep = "")
	tt=fdataobj$argvals
	rtt=fdataobj$rangval
    n <- nrow(fdataobj)
    m <- ncol(fdataobj)
    m2 <- ncol(fdataori)
    n2 <- nrow(fdataori)
			}
    else {
        stop("no fdata class object on input")
    }
    if (is.null(n) && is.null(m)) 
        stop("ERROR IN THE DATA DIMENSIONS")
    if (is.null(m) && is.null(m2)) 
        stop("ERROR IN THE DATA DIMENSIONS")
    mdist=matrix(NA,ncol=n2,nrow=n2)
    for (i in 1:(n2-1)){for (j in (i+1):n2){
      mdist[i,j]<-mdist[j,i]<-norm.fdata(fdataori[i]-fdataori[j])
    }}
    if (is.null(h))   {
	  h<-0.15
	  hq2=quantile(mdist,probs=h,na.rm=TRUE)
	  #print("es nulo h")  
	}
	else {
	  #cat("no es nulo h ",h,"\n")    
	  if (is.numeric(h))    hq2<-h  
	  else hq2=quantile(mdist,probs=as.numeric(h),na.rm=TRUE)
	}	
	kern=function(x,y,h=hq2){exp(-norm.fdata(x-y)^2/h^2)}
	K02=rep(1,nrow(fdataobj))
	K01=rep(1,nrow(fdataori))
	M1=matrix(NA,nrow=n2,ncol=n2)
	M2=matrix(NA,nrow=n,ncol=n2)
	M=array(NA,dim=c(n,n2,n2))
	for (i in 1:n2){for (j in i:n2){
	if (i==j) M1[i,i]=K01[i] else M1[i,j]<-M1[j,i]<-kern(fdataori[i],fdataori[j],h=hq2)
	}}
	same.dim <- FALSE
	if (n==n2 & m==m2){ same.dim <- TRUE}
	if (same.dim)
	  if (all(fdataobj$data==fdataori$data)) {
	    M2=M1
	  } else {sam.dim <- FALSE}
	if (!same.dim){	
  	for (i in 1:n){for (j in 1:n2){
  	if (all(fdataobj[i]$data == fdataori[j]$data)) M2[i,j]=K02[i] else M2[i,j]=kern(fdataobj[i],fdataori[j],h=hq2)
  	}}
	}  
	for (i in 1:n){for (j in 1:n2){for (k in 1:n2){
	M[i,j,k]<-(K02[i]+M1[j,k]-M2[i,j]-M2[i,k])/(sqrt(K02[i]+K01[j]-2*M2[i,j])*sqrt(K02[i]+K01[k]-2*M2[i,k]))
	}}}
	l=which(!is.finite(M))
	M[l]=NA
	dep=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n
    if (scale) {
	MO=array(NA,dim=c(n2,n2,n2))
	for (i in 1:n2){for (j in 1:n2){for (k in 1:n2){
	MO[i,j,k]<-(K01[i]+M1[j,k]-M1[i,j]-M1[i,k])/(sqrt(K01[i]+K01[j]-2*M1[i,j])*sqrt(K01[i]+K01[k]-2*M1[i,k]))
	}}}
	l=which(!is.finite(MO))
	MO[l]=NA
	dep2=1-sqrt(apply(MO,1,sum,na.rm=TRUE))/n2
        mn <- min(dep2, na.rm = TRUE)
        mx <- max(dep2, na.rm = TRUE)
        dep = as.vector(dep/mx)
    }
    if (nullans) {
        ans1 <- rep(NA, len = m0)
        ans1[-nas] <- dep
        dep <- ans1
    }
	names(dep)=nms
    k = which.max(dep)
    med = fdataobj[k]
    nl = length(trim)
	lista=vector("list",nl)
    tr <- paste("KFSD.tr", round(trim * 100,2), "%", sep = "")
    if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
    mtrim = fdata(matrix(NA, nrow = nl, ncol = m),tt,rtt,names)
    for (j in 1:length(trim)) {
        lista[[j]] = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
        if (length(lista[[j]])==1) {
          mtrim$data[j,]<-fdataobj[lista[[j]]]$data
          if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
        }
        else mtrim$data[j,]=func.mean(fdataobj[lista[[j]]])$data       
    }
    rownames(med$data) <- "KFSD.med"
    rownames(mtrim$data) <- tr
    out <- list(median = med, lmed = k, mtrim = mtrim, ltrim = if (nl==1) unlist(lista) else lista, 
        dep = dep, h = h, hq=hq2,name="KFSD")
    if (scale) out$dscale = mx
	class(out)="depth"
	out$trim <- trim
	out$name <- "KFSD"
	out$fdataobj=fdataobj
	out$fdataori=fdataori
    if (draw) {
	    plot.depth(out)
    }
    return(invisible(out))	
}
#####################
depth.FSD=function (fdataobj, fdataori = fdataobj, 
                    trim = 0.25, scale = FALSE, draw = FALSE){
    if (is.fdata(fdataobj)) {
        if (is.null(rownames(fdataobj$data))) 
            rownames(fdataobj$data) <- 1:nrow(fdataobj$data)
        m0=nrow(fdataobj)
        nms <- rownames(fdataobj$data)
        fdataobj <- na.omit.fdata(fdataobj)
        fdataori <- na.omit.fdata(fdataori)
        nas <- na.action(fdataobj)
		n <- nrow(fdataobj)
		m <- ncol(fdataobj)
		m2 <- ncol(fdataori)
		n2 <- nrow(fdataori)
        nullans <- !is.null(nas)
        names1 <- names2 <- names <- fdataobj[["names"]]
        names1$main <- "depth.FSD median"
        names2$main <- paste("depth.FSD trim ", trim * 100, 
            "%", sep = "")
        tt = fdataobj[["argvals"]]
        rtt <- fdataobj[["rangeval"]]
    }
    else {
        stop("no fdata class object")
    }
    if (is.null(n) && is.null(m)) 
        stop("ERROR IN THE DATA DIMENSIONS")
    if (is.null(m) && is.null(m2)) 
        stop("ERROR IN THE DATA DIMENSIONS")
	kern=function(x,y){inprod.fdata(x,y)}

	K02=norm.fdata(fdataobj)^2
	K01=norm.fdata(fdataori)^2
	M1=matrix(NA,nrow=n2,ncol=n2)
	M2=matrix(NA,nrow=n,ncol=n2)
	M=array(NA,dim=c(n,n2,n2))
	if (scale) MO=array(NA,dim=c(n2,n2,n2))
	for (i in 1:n2){for (j in i:n2){
	if (i==j) M1[i,i]=K01[i] else M1[i,j]<-M1[j,i]<-kern(fdataori[i],fdataori[j])
	}}
  same.dim <- FALSE
  if (n==n2 & m==m2){ same.dim <- TRUE}
    if (same.dim)
      if (all(fdataobj$data==fdataori$data)) {
    	 M2=M1
	    } else {sam.dim <- FALSE}
  if (!same.dim){
  	for (i in 1:n){for (j in 1:n2){
  	if (all(fdataobj[i]$data == fdataori[j]$data)) M2[i,j]=K02[i] else M2[i,j]=kern(fdataobj[i],fdataori[j])
  }}
	}
	for (i in 1:n){for (j in 1:n2){for (k in 1:n2){
	if (all(fdataobj[i]$data == fdataori[j]$data) | all(fdataobj[i]$data == fdataori[k]$data)) M[i,j,k]=NA else M[i,j,k]<-(K02[i]+M1[j,k]-M2[i,j]-M2[i,k])/(sqrt(K02[i]+K01[j]-2*M2[i,j])*sqrt(K02[i]+K01[k]-2*M2[i,k]))
	}}}
	l=which(!is.finite(M))
	M[l]=NA
	dep=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n
    if (scale) {
	MO=array(NA,dim=c(n2,n2,n2))
	for (i in 1:n2){for (j in 1:n2){for (k in 1:n2){
		if (all(fdataori[i]$data == fdataori[j]$data) | all(fdataori[i]$data == fdataori[k]$data)) MO[i,j,k]=NA else MO[i,j,k]<-(K01[i]+M1[j,k]-M1[i,j]-M1[i,k])/(sqrt(K01[i]+K01[j]-2*M1[i,j])*sqrt(K01[i]+K01[k]-2*M1[i,k]))
	}}}
	l=which(!is.finite(MO))
	MO[l]=NA
	dep2=1-sqrt(apply(MO,1,sum,na.rm=TRUE))/n2
	mn <- min(dep2, na.rm = TRUE)
    mx <- max(dep2, na.rm = TRUE)
    dep = as.vector(dep/mx)
    }
    if (nullans) {
        ans1 <- rep(NA, len = m0)
        ans1[-nas] <- dep
        dep <- ans1
    }
	names(dep)=nms
    k = which.max(dep)
    med = fdataobj[k]
    nl = length(trim)
	lista=vector("list",nl)
	tr <- paste("FSD.tr", round(trim * 100,2), "%", sep = "")
    if (nl>1) names(lista)=paste0("tr",round(trim*100,2))	
    mtrim = fdata(matrix(NA, nrow = nl, ncol = m),tt,rtt,names2)
    for (j in 1:length(trim)) {
        lista[[j]] = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
        if (length(lista[[j]])==1) {
          mtrim$data[j,]<-fdataobj[lista[[j]]]$data
          if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
        }
        else mtrim$data[j,]=func.mean(fdataobj[lista[[j]]])$data       
    }
  
    rownames(med$data) <- "FSD.med"
    rownames(mtrim$data) <- tr
    out <- list(median = med, lmed = k, mtrim = mtrim, ltrim = if (nl==1) unlist(lista) else lista, 
        dep = dep,name="FSD")
    if (scale) 
        out$dscale = mx
	 class(out)="depth"
	 out$trim <- trim
	 out$name <- "FSD"
	 out$fdataobj=fdataobj
	 out$fdataori=fdataori
    if (draw) {
	plot.depth(out)
    }
     return(invisible(out))
	}

####################
mdepth.KFSD=function (dataobj, dataori = dataobj, trim = 0.25,
                     h=NULL,scale = FALSE, draw = FALSE){
	if (is.matrix(dataobj) & is.matrix(dataori)){
		m0=nrow(dataobj)
        rownames(dataobj) <- 1:nrow(dataobj)
        nms <- rownames(dataobj)
        dataobj <- na.omit(dataobj)
        dataori <- na.omit(dataori)
		m2 <- ncol(dataori)
		n2 <- nrow(dataori)	
		n <- nrow(dataobj)
		m <- ncol(dataobj)		
	if (m2!=m) stop ("Error in dimensions of dataobj and dataori")
        nas <- na.action(dataobj)
        nullans <- !is.null(nas)
    }
    else {
        stop("no matrix object")
    }
 
    if (is.null(n) && is.null(m)) 
        stop("ERROR IN THE DATA DIMENSIONS")
    if (is.null(m) && is.null(m2)) 
        stop("ERROR IN THE DATA DIMENSIONS")
    mdist=matrix(NA,ncol=n2,nrow=n2)
    for (i in 1:(n2-1)){for (j in (i+1):n2){
      mdist[i,j]<-mdist[j,i]<-sqrt(sum((dataori[i,]-dataori[j,])^2))
    }}
    if (is.null(h))   {
	  h<-0.15
	  hq2=quantile(mdist,probs=h,na.rm=TRUE)
	  #print("es nulo h")  
	}
	else {
	  #cat("no es nulo h ",h,"\n")    
	  if (is.numeric(h))    hq2<-h  
	  else hq2=quantile(mdist,probs=as.numeric(h),na.rm=TRUE)
	}	
	kern=function(x,y,h=hq2){exp(-sum((x-y)^2)/h^2)}
	K02=rep(1,nrow(dataobj))
	K01=rep(1,nrow(dataori))
	M1=matrix(NA,nrow=n2,ncol=n2)
	M2=matrix(NA,nrow=n,ncol=n2)
	M=array(NA,dim=c(n,n2,n2))
	for (i in 1:n2){for (j in i:n2){
	if (i==j) M1[i,i]=K01[i] else M1[i,j]<-M1[j,i]<-kern(dataori[i,],dataori[j,],h=hq2)
	}}
	same.dim <- FALSE
	if (n==n2 & m==m2){ same.dim <- TRUE}
	if (same.dim)
	  if (all(dataobj==dataori)) {
	    M2=M1
	  } else {sam.dim <- FALSE}
	if (!same.dim){	
  	for (i in 1:n){for (j in 1:n2){
  	if (all(dataobj[i,] == dataori[j,])) M2[i,j]=K02[i] else M2[i,j]=kern(dataobj[i,],dataori[j,],h=hq2)
  	}}}  
	for (i in 1:n){for (j in 1:n2){for (k in 1:n2){	
	M[i,j,k]<-(K02[i]+M1[j,k]-M2[i,j]-M2[i,k])/(sqrt(K02[i]+K01[j]-2*M2[i,j])*sqrt(K02[i]+K01[k]-2*M2[i,k]))
	}}}
	l=which(!is.finite(M))
	M[l]=NA
	dep=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n
    if (scale) {
	MO=array(NA,dim=c(n2,n2,n2))
	for (i in 1:n2){for (j in 1:n2){for (k in 1:n2){	
	MO[i,j,k]<-(K02[i]+M1[j,k]-M1[i,j]-M1[i,k])/(sqrt(K02[i]+K01[j]-2*M1[i,j])*sqrt(K02[i]+K01[k]-2*M1[i,k]))
	}}}
	l=which(!is.finite(MO))
	MO[l]=NA
	dep2=1-sqrt(apply(MO,1,sum,na.rm=TRUE))/n2
	mn <- min(dep2, na.rm = TRUE)
    mx <- max(dep2, na.rm = TRUE)
    dep = as.vector(dep/mx)
    }
   if (nullans) {
        ans1 <- rep(NA, len = m0)
        ans1[-nas] <- dep
        dep <- ans1
    }
	names(dep)=nms
    k = which.max(dep)
    med = matrix(dataobj[k,],nrow=1)
    nl = length(trim)
	lista=vector("list",nl)
    tr <- paste("KFSD.tr", round(trim * 100,2), "%", sep = "")
	if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
    mtrim = matrix(NA, nrow = nl, ncol = m)
    for (j in 1:length(trim)) {
        lista[[j]] = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
        if (length(lista[[j]])==1) {
          mtrim[j,]<-dataobj[lista[[j]],]
          if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
        }
        else mtrim[j,]=apply(dataobj[lista[[j]],],2,mean)       
    }

    rownames(med) <- "KFSD.med"
    rownames(mtrim) <- tr
    out <- list(median = med, lmed = k, mtrim = mtrim, ltrim = if (nl==1) unlist(lista) else lista, 
        dep = dep, h = h, hq=hq2,name="KFSD")
    if (scale) 
        out$dscale = mx
	class(out)="mdepth"
    if (draw) {
		plot.mdepth(dataobj,dataori,out)
    }
    return(invisible(out))	
}

mdepth.FSD=function (dataobj, dataori = dataobj, 
                    trim = 0.25, scale = FALSE, draw = FALSE){
    if (is.matrix(dataobj) & is.matrix(dataori)) {
        if (is.null(rownames(dataobj))) 
            rownames(dataobj) <- 1:nrow(dataobj)
        nms <- rownames(dataobj)
        m0 <- nrow(dataobj)
        dataobj <- na.omit(dataobj)
        dataori <- na.omit(dataori)
        nas <- na.action(dataobj)
        nullans <- !is.null(nas)
		n <- nrow(dataobj)
		m <- ncol(dataobj)
		m2 <- ncol(dataori)
		n2 <- nrow(dataori)
		if (m2!=m) stop("Error in data dimensions")
		}
    else {
        stop("no matrix object")
    }
    if (is.null(n) && is.null(m)) 
        stop("ERROR IN THE DATA DIMENSIONS")
    if (is.null(m) && is.null(m2)) 
        stop("ERROR IN THE DATA DIMENSIONS")
	kern=function(x,y){sum(x*y)}
	K02=diag(dataobj%*%t(dataobj))
	K01=diag(dataori%*%t(dataori))
	M1=matrix(NA,nrow=n2,ncol=n2)
	M2=matrix(NA,nrow=n,ncol=n2)
	M=array(NA,dim=c(n,n2,n2))
	if (scale) MO=array(NA,dim=c(n2,n2,n2))
	for (i in 1:n2){for (j in i:n2){
	if (i==j) M1[i,i]=K01[i] else M1[i,j]<-M1[j,i]<-kern(dataori[i,],dataori[j,])
	}}
  same.dim <- FALSE
  if (n==n2 & m==m2){ same.dim <- TRUE}
    if (same.dim)
      if (all(dataobj==dataori)) {
    	 M2=M1
	    } else {sam.dim <- FALSE}
  if (!same.dim){
  	for (i in 1:n){for (j in 1:n2){
  	if (all(dataobj[i,] == dataori[j,])) M2[i,j]=K02[i] else M2[i,j]=kern(dataobj[i,],dataori[j,])
  }}
	}
	for (i in 1:n){for (j in 1:n2){for (k in 1:n2){
	if (all(dataobj[i,] == dataori[j,]) | all(dataobj[i,] == dataori[k,])) M[i,j,k]=NA else M[i,j,k]<-(K02[i]+M1[j,k]-M2[i,j]-M2[i,k])/(sqrt(K02[i]+K01[j]-2*M2[i,j])*sqrt(K02[i]+K01[k]-2*M2[i,k]))
	}}}
	l=which(!is.finite(M))
	M[l]=NA
	dep=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n
    if (scale) {
		for (i in 1:n2){for (j in 1:n2){for (k in 1:n2){
		if (all(dataori[i,] == dataori[j,]) | all(dataori[i,] == dataori[k,])) MO[i,j,k]=NA else MO[i,j,k]<-(K01[i]+M1[j,k]-M1[i,j]-M1[i,k])/(sqrt(K01[i]+K01[j]-2*M1[i,j])*sqrt(K01[i]+K01[k]-2*M1[i,k]))
	}}}
	l=which(!is.finite(MO))
		MO[l]=NA
	dep2=1-sqrt(apply(MO,1,sum,na.rm=TRUE))/n2
		mn <- min(dep2, na.rm = TRUE)
        mx <- max(dep2, na.rm = TRUE)
        dep = as.vector(dep/mx)
    }
   if (nullans) {
        ans1 <- rep(NA, len = m0)
        ans1[-nas] <- dep
        dep <- ans1
    }
	names(dep)=nms
    k = which.max(dep)
    med = matrix(dataobj[k,],nrow=1)
    nl = length(trim)
	lista=vector("list",nl)
    tr <- paste("FSD.tr", round(trim * 100,2), "%", sep = "")
	if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
    mtrim = matrix(NA, nrow = nl, ncol = m)
    for (j in 1:length(trim)) {
        lista[[j]] = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
        if (length(lista[[j]])==1) {
          mtrim[j,]<-dataobj[lista[[j]],]
          if (draw) {draw=FALSE;warning("Too data in mtrim. The plot is not drawn")}
        }
        else mtrim[j,]=apply(dataobj[lista[[j]],],2,mean)       
    }

    rownames(med) <- "FSD.med"
    rownames(mtrim) <- tr
    out <- list(median = med, lmed = k, mtrim = mtrim, ltrim = if (nl==1) unlist(lista) else lista, 
        dep = dep,name="FSD")
    if (scale) 
        out$dscale = mx
	class(out)="mdepth"
	
		if (draw) {
		plot.mdepth(out)
			}
    return(invisible(out))
	}

