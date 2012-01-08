 depth.RT  <-function (fdataobj,fdataori=fdataobj, trim = 0.25, nproj = 10, proj = 1, xeps = 1e-07, 
    draw = FALSE, ...) 
{
    if (!is.fdata(fdataobj)) fdataobj = fdata(fdataobj)
    if (!is.fdata(fdataori)) fdataobj=fdata(fdataori)    
    nas <- apply(fdataobj$data, 1, count.na)
    if (any(nas)) {
        fdataobj$data <- fdataobj$data[!nas, ]
        cat("Warning: ", sum(nas), " curves with NA are not used in the calculations \n")
    }
    data  <- fdataobj[["data"]]
    data2 <- fdataori[["data"]]    
    n <- nrow(data)
    m <- ncol(data)
    m2<-ncol(data2)
    n2<-nrow(data2)
    if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS OF fdataobj")
    if (is.null(n2) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS fdataori")
    t = fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names1 <- names2 <- names <- fdataobj[["names"]]
    names1$main <- "depth.RP median"
    names2$main <- paste("RP trim", trim * 100, "%", sep = "")
   
    if (is.fdata(proj)) {
        nproj <- nrow(proj)
        if (fdataobj$argvals != proj$argvals || ncol(fdataobj) != 
            ncol(proj)) 
            stop("Error en proj dimension")
        z <- proj
    }
    else {
        z <- rproc2fdata(nproj, t, sigma = proj, norm = TRUE,...)
    }
    Fn <- list()
    Prod=data%*%t(z$data)
    Prod2=data2%*%t(z$data)    
    dep = array(NA, dim = c(nproj, n))
    dep2 = array(NA, dim = c(nproj, n2))    
    
 #   dep2=rep(0.0,n) 
    for (j in 1:nproj) {
        Fn[[j]] = ecdf(Prod2[,j])  
        dep[j, ] = pmin(Fn[[j]](Prod[,j]) ,(1 - Fn[[j]](Prod[,j]-xeps)))       
        dep2[j, ] = pmin(Fn[[j]](Prod2[,j]) ,(1 - Fn[[j]](Prod2[,j]-xeps)))               
#        dep2 =dep2+ pmin(Fn[[j]](Prod[,j]) ,(1 - Fn[[j]](Prod[,j]))) #idem que dep
    }
#    print(dep)
#    dep = apply(dep, 2, min)
    dep = colMeans(dep)
    dep2 = colMeans(dep2)    
#    print(cbind(dep,dep2/nproj,dep3))
    k = which.max(dep2)
    med = data2[k, ]
    nl = length(trim)
    mtrim = matrix(NA, nrow = nl, ncol = m)
    for (j in 1:length(trim)) {
        lista = which(dep2 >= quantile(dep2, probs = trim[j], na.rm = TRUE))
        if (length(lista)==1) {
          mtrim[j, ]<-data2[lista,]
          if (draw) {draw=FALSE;warning("The plot is not shown")}
        }
        else mtrim[j, ]=apply(data2[lista,],2,mean,na.rm=TRUE)       
    }
    tr <- paste("RP.tr", trim * 100, "%", sep = "")
    med <- fdata(med, t, rtt, names1)
    mtrim <- fdata(mtrim, t, rtt, names2)
    rownames(med$data) <- "RP.med"
    rownames(mtrim$data) <- tr
    if (draw) {
        ans <- dep2
        ind1 <- !is.na(ans)
        cgray = 1 - (ans - min(ans, na.rm = TRUE))/(max(ans, 
            na.rm = TRUE) - min(ans, na.rm = TRUE))
        plot(fdataori[ind1, ], col = gray(cgray[ind1]), main = "RP Depth")
        lines(mtrim, lwd = 2, col = "yellow")
        lines(med, col = "red", lwd = 2)
        legend("topleft", legend = c(tr, "Median"), lwd = 2,box.col=0, 
            col = c("yellow", "red"))
    }
    return(invisible(list(median = med, lmed = k, mtrim = mtrim, 
        ltrim = lista, dep = dep,dep.ori=dep, proj = z, Fn = Fn)))
}




