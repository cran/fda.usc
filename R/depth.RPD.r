depth.RPD<-function (fdataobj,fdataori=fdataobj, nproj = 20, proj=1,deriv = c(0, 1), trim = 0.25,
    dfunc2 = mdepth.LD, method = "fmm", draw = FALSE, ...)
{
    if (!is.fdata(fdataobj))         fdataobj = fdata(fdataobj)
    if (!is.fdata(fdataori)) fdataobj=fdata(fdataori)    
	if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
 nms<-rownames(fdataobj$data)
 m0<-nrow(fdataobj)
 fdataobj2<- fdataobj
 fdataobj<-na.omit.fdata(fdataobj)
 fdataori<-na.omit.fdata(fdataori) 
 nas<-na.action(fdataobj)
 nullans<-!is.null(nas) 
    names1 <- names2 <- fdataobj[["names"]]
    names1$main <- "depth.RPD median"
    names2$main <- paste("RPD trim ", trim * 100, "%", sep = "")
    n <- nrow(fdataobj)
    m <- ncol(fdataobj)
    m2<-ncol(fdataori)
    n2<-nrow(fdataori)    
    if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
    if (is.null(n) || is.null(m))       stop("Input must be a matrix")
    tt = fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
	newfunc=vector("list",length(deriv))
	newfunc2=vector("list",length(deriv))
    for (ider in 1:length(deriv)) {
        if (deriv[ider] == 0) {   
		  newfunc[[ider]]=fdataobj
		  newfunc2[[ider]]=fdataori
          } 
        else {
          newfunc[[ider]] = fdata.deriv(fdataobj, nderiv = deriv[ider], method = method, ...)
          newfunc2[[ider]] = fdata.deriv(fdataori, nderiv = deriv[ider], method = method, ...) 
        }
    }
    dep = rep(0, n)
    dep2 = rep(0, n2)    
    vproject = matrix(0, nrow = n, ncol = length(deriv))
    vproject2 = matrix(0, nrow = n2, ncol = length(deriv))    
    if (is.fdata(proj)) {
     if (fdataobj$argvals!=proj$argvals || m!=ncol(proj)) stop("Error in proj dimension")
     z<-proj
     nproj<-nrow(z)
    }
  else {	 z<-rproc2fdata(nproj,tt,sigma=proj,norm=TRUE,...)	}
  pb = txtProgressBar(min = 0, max = nproj, style = 3)
  for (j in 1:nproj) {
        setTxtProgressBar(pb, j - 0.5)
        for (ider in 1:length(deriv)) {
            vproject[, ider] = inprod.fdata(newfunc[[ider]],z[j])
            vproject2[, ider] = inprod.fdata(newfunc2[[ider]],z[j]) 
        }
        par.dfunc = list()
        par.dfunc$x <- vproject
        par.dfunc$xx <- vproject2
#        par.dfunc$trim <- trim              
        par.dfunc$scale<-TRUE
        resul = do.call(dfunc2, par.dfunc)
        dep = dep + resul$dep
        setTxtProgressBar(pb, j)
    }                                                               
    close(pb)
	names(dep)<-nms       
    dep = dep/nproj
if (nullans) {
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-dep
        dep<-ans1      
        }
    k = which.max(dep)    
    med = fdataobj[k]
   nl=length(trim)
   lista=vector("list",nl)
   tr<-paste("RPD.tr",round(trim*100,2),"\u0025",sep="")
   if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
   mtrim=fdata(matrix(NA,ncol=length(tt),nrow=length(trim)),tt,rtt,names2)
   for (j in 1:nl){
   lista[[j]]=which(dep>=quantile(dep,probs=trim[j],na.rm=TRUE))
   	if (length(lista[[j]])==1) {
			mtrim$data[j,]<-fdataobj2[lista[[j]]]$data
	if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
	}
	else   mtrim$data[j,]=apply(fdataobj2$data[lista[[j]],,drop=FALSE],2,mean,na.rm=TRUE)
	}
#else mtrim=data[lista,,drop=FALSE]
#   med<-fdata(med,tt,rtt,names1)
#   mtrim<-fdata(mtrim,tt,rtt,names2)
   rownames(med$data)<-"RPD.med"
   rownames(mtrim$data)<-tr
   out=list(median = med, lmed = k, mtrim = mtrim,
        ltrim = if (nl==1) unlist(lista) else lista, dep = dep,deriv=deriv,proj = z,name="RPD")
   out$trim <- trim
   out$name <- "RPD"
   out$fdataobj=fdataobj
   out$fdataori=fdataori
	class(out)="depth"
    if (draw) {
	  plot.depth(out)
		}
    return(invisible(out))
}  
