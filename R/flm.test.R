
##############################################
##          Functional Linear Model         ##
##             Goodness-of-fit              ##
##############################################

##############################################
## File created by Eduardo García-Portugués ##
## using code from library fda.usc 			##
##############################################

# Generation of centred random variables with unit variance for the wild bootstrap
rber.gold=function(n){
	res=rbinom(n=n,size=1,prob=(5+sqrt(5))/10)
	res=ifelse(res==1,(1-sqrt(5))/2,(1+sqrt(5))/2)
	return(res)
}

# Aijr function to call Fortran code
Aijr=function(X,inpr){
	
	# Check if the Fortran function is correctly loaded
	if(!is.loaded("aijr")) stop("Fortran function aijr is not loaded")
	
	# Compute the inproduct
	if(missing(inpr)){
		if(is.fd(X)){
			inpr=t(X$coefs)%*%inprod(X$basis,X$basis)%*%X$coefs
		}else if(is.fdata(X)){
			inpr=inprod.fdata(X,X)
		}else{
			stop("X is not a fd or fdata object")
		}
	}
	
	# Number of functions in X
	n=dim(inpr)[1]
	
	# Check for repeated interior products (will cause error in the aijr subroutine)
	if(anyDuplicated(diag(inpr))){
		diag(inpr)=diag(inpr)*(1+rexp(n,rate=1000))
	}
	
	# Call Fortran function
	res=.Fortran("aijr",n=as.integer(n),inprod=matrix(as.double(inpr),nrow=n,ncol=n),Aijr_vec=numeric((n^3+n^2)/2))
	
	# Return result
	return(res$Aijr_vec)
	
}

# PCvM statistic
PCvM.statistic=function(X,residuals,p,Aijr0){
	
	# Check if the Fortran function is correctly loaded
	if(!is.loaded("pcvm_statistic")) stop("Fortran function pcvm_statistic is not loaded")
	
	# Obtain the sample size
	n=length(residuals)
	
	# Check if the lengths of e and X are the same
	if(is.fd(X)){
		if(n!=dim(X$coefs)[2]) stop("Incompatible lenghts of X.fd and residuals")
	}else if(is.fdata(X)){
		if(n!=dim(X$data)[1]) stop("Incompatible lenghts of X.fdata and residuals")
	}
	
	
	# Compute Aijr if missing
	if(missing(Aijr0)) Aijr0=Aijr(X)
	
	# Call Fortran function
	res=n^(-2)*pi^((p/2)-1)/gamma(p/2+1)*.Fortran("pcvm_statistic",n=as.integer(n),Aijr_vec=Aijr0,residuals=residuals,statistic=0)$statistic
	
	return(res)

}

# PCvM test for the composite hypothesis with bootstrap calibration
flm.test=function(X.fdata,Y,beta0.fdata=NULL,B=5000,est.method="pls",p=NULL,type.basis="bspline",show.prog=TRUE,plot.it=TRUE,B.plot=100,G=200,...){
	
	# Center the data first
	X.fdata=fdata.cen(X.fdata)$Xcen
	Y=Y-mean(Y)
	
	# Number of functions
	n=dim(X.fdata)[1]
	
	# Bootstrap residuals
	if(B.plot>B) stop("B.plot must be less or equal than B")
	e.hat.star=matrix(ncol=n,nrow=B)
	boot.beta.est=list()
	
	## COMPOSITE HYPOTHESIS ##
	if(show.prog) cat("Computing estimation of beta... ")
	if(is.null(beta0.fdata)){

		## 1. Optimal estimation of beta and the basis order ##
		
		if(est.method=="pc"){
			
			if(is.null(p)){

				# Method
				meth="PCvM test for the functional linear model using optimal PC basis representation"

				# Choose the number of basis elements: SIC is probably the best criteria
				mod.pc=fregre.pc.cv(fdataobj=X.fdata,y=Y,kmax=10,criteria="SIC")
				p.opt=length(mod.pc$pc.opt)
				ord.opt=mod.pc$pc.opt
				
				# PC components to be passed to the bootstrap
				pc.comp=mod.pc$fregre.pc$fdata.comp # pc.comp=mod.pc$fregre.pc$pc
				pc.comp$l=mod.pc$pc.opt
				
				# Express X.fdata and beta.est in the PC basis
				basis.pc=mod.pc$fregre.pc$fdata.comp$rotation
				if(length(pc.comp$l)!=1){
					X.est=fdata(mdata=mod.pc$fregre.pc$fdata.comp$x[,mod.pc$fregre.pc$l]%*%mod.pc$fregre.pc$fdata.comp$rotation$data[mod.pc$fregre.pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval) # X.est=fdata(mdata=mod.pc$fregre.pc$pc$x[,mod.pc$fregre.pc$l]%*%mod.pc$fregre.pc$pc$rotation$data[mod.pc$fregre.pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.pc$fregre.pc$fdata.comp$x[,mod.pc$fregre.pc$l]%*%t(mod.pc$fregre.pc$fdata.comp$rotation$data[mod.pc$fregre.pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval) # X.est=fdata(mdata=mod.pc$fregre.pc$pc$x[,mod.pc$fregre.pc$l]%*%t(mod.pc$fregre.pc$pc$rotation$data[mod.pc$fregre.pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.pc$fregre.pc$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.pc$fregre.pc$residuals
				
			}else{
			
				# Method
				meth=paste("PCvM test for the functional linear model using a representation in a PC basis of ",p,"elements") 
				
				# Estimation of beta on the given fixed basis
				mod.fixed=fregre.pc(fdataobj=X.fdata,y=Y,l=1:p)
				p.opt=p
				ord.opt=mod.fixed$l

				# PC components to be passed to the bootstrap
				pc.comp=mod.fixed$pc
				pc.comp$l=mod.fixed$l
				
				# Express X.fdata and beta.est in the basis
				if(p!=1){
					X.est=fdata(mdata=mod.fixed$fdata.comp$x[,mod.fixed$l]%*%mod.fixed$fdata.comp$rotation$data[mod.fixed$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.fixed$fdata.comp$x[,mod.fixed$l]%*%t(mod.fixed$fdata.comp$rotation$data[mod.fixed$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.fixed$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.fixed$residuals
			
			}
			
		}else if(est.method=="pls"){
			
			if(is.null(p)){

				# Method
				meth="PCvM test for the functional linear model using optimal PLS basis representation"
			
				# Choose the number of the basis: SIC criteria overfit the data, we choose CV althought is slower
				mod.pls=fregre.plsr.cv(fdataobj=X.fdata,y=Y,kmax=10,criteria="CV") 
				p.opt=length(mod.pls$pls.opt)
				ord.opt=mod.pls$pls.opt
				
				# PLS components to be passed to the bootstrap
				pls.comp=mod.pls$fregre.pls$fdata.comp
				pls.comp$l=mod.pls$pls.opt
						
				# Express X.fdata and beta.est in the PLS basis
				basis.pls=mod.pls$fregre.pls$fdata.comp$rotation
				if(length(pls.comp$l)!=1){
					X.est=fdata(mdata=mod.pls$fregre.pls$fdata.comp$x[,mod.pls$fregre.pls$l]%*%mod.pls$fregre.pls$fdata.comp$rotation$data[mod.pls$fregre.pls$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.pls$fregre.pls$fdata.comp$x[,mod.pls$fregre.pls$l]%*%t(mod.pls$fregre.pls$fdata.comp$rotation$data[mod.pls$fregre.pls$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.pls$fregre.pls$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.pls$fregre.pls$residuals
				
			}else{
			
				# Method
				meth=paste("PCvM test for the functional linear model using a representation in a PLS basis of ",p,"elements") 
			
				# Estimation of beta on the given fixed basis
				mod.fixed=fregre.pc(fdataobj=X.fdata,y=Y,l=1:p)
				p.opt=p
				ord.opt=mod.fixed$l

				# PLS components to be passed to the bootstrap
				pls.comp=mod.fixed$fdata.comp
				pls.comp$l=mod.fixed$l
				
				# Express X.fdata and beta.est in the basis
				if(p!=1){
					X.est=fdata(mdata=mod.fixed$fdata.comp$x[,mod.fixed$l]%*%mod.fixed$fdata.comp$rotation$data[mod.fixed$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.fixed$fdata.comp$x[,mod.fixed$l]%*%t(mod.fixed$fdata.comp$rotation$data[mod.fixed$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.fixed$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.fixed$residuals
			
			}
			
		}else if(est.method=="basis"){
			
			if(is.null(p)){
			
				# Method
				meth=paste("PCvM test for the functional linear model using optimal",type.basis,"basis representation")
			
				# Choose the number of the bspline basis with GCV.S
				mod.basis=fregre.basis.cv.mod(x=X.fdata,y=Y,p=seq(5,30,by=1),type.basis=type.basis,type.CV=GCV.S,...)
				p.opt=mod.basis$basis.x.opt$nbasis
				ord.opt=1:p.opt
				
				# Express X.fdata and beta.est in the optimal basis
				basis.opt=mod.basis$basis.x.opt
				X.est=mod.basis$x.fd
				beta.est=mod.basis$beta.est
				norm.beta.est=norm.fd(beta.est)
				
				# Compute the residuals
				e=mod.basis$residuals
			
			}else{
				
				# Method
				meth=paste("PCvM test for the functional linear model using a representation in a",type.basis,"basis of ",p,"elements") 

				# Estimation of beta on the given fixed basis
				basis.opt=do.call(what=paste("create.",type.basis,".basis",sep=""),args=list(rangeval=X.fdata$rangeval,nbasis=p,...))
				mod.fixed=fregre.basis(fdataobj=X.fdata,y=Y,basis.x=basis.opt,basis.b=basis.opt)
				p.opt=p
				ord.opt=1:p.opt

				# Express X.fdata and beta.est in the basis
				X.est=mod.fixed$x.fd
				beta.est=mod.fixed$beta.est
				norm.beta.est=norm.fd(beta.est)
				
				# Compute the residuals
				e=mod.fixed$residuals
			
			}
			
		}else{
		
			stop(paste("Estimation method",est.method,"not implemented."))
			
		}
	
	## SIMPLE HYPOTHESIS ##
	}else{
	
		## 1. Optimal representation of X and beta0 ##
		
		# Choose the number of basis elements
		if(type.basis!="pc" & type.basis!="pls"){
			
			# Basis method: select the number of elements of the basis is it is not given
			if(is.null(p)){
				
				# Method
				meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in an optimal ",type.basis,"basis") 
				
				p.opt=min.basis.mod(X.fdata,type.basis=type.basis,numbasis=seq(31,71,by=2))$numbasis.opt
				if(p.opt==31){
					cat("Readapting interval...\n")
					p.opt=min.basis.mod(X.fdata,type.basis=type.basis,numbasis=seq(5,31,by=2))$numbasis.opt
				}else if(p.opt==71){
					cat("Readapting interval...\n")
					p.opt=min.basis.mod(X.fdata,type.basis=type.basis,numbasis=seq(71,101,by=2))$numbasis.opt
				}
				
				ord.opt=1:p.opt
				
			}else{
				
				# Method
				meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in a ",type.basis," basis of ",p,"elements") 
				
				p.opt=p
				ord.opt=1:p.opt
			}
			
			# Express X.fdata in the basis
			X.est=fdata2fd(X.fdata,type.basis=type.basis,nbasis=p)
			beta.est=fdata2fd(beta0.fdata,type.basis=type.basis,nbasis=p)
			
			# Compute the residuals
			e=Y-inprod(X.est,beta.est)

		}else if (type.basis=="pc"){
			
			if(is.null(p)) stop("Simple hypothesis with type.basis=\"pc\" need the number of components p") 
			# Method
			meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in a PC basis of ",p,"elements") 
			
			# Express X.fdata in a PC basis
			fd2pc=fdata2pc(X.fdata,ncomp=p)
			if(length(fd2pc$l)!=1){
				X.est=fdata(mdata=fd2pc$x[,fd2pc$l]%*%fd2pc$rotation$data[fd2pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}else{
				X.est=fdata(mdata=fd2pc$x[,fd2pc$l]%*%t(fd2pc$rotation$data[fd2pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}
			beta.est=beta0.fdata
			p.opt=p
			ord.opt=1:p.opt
			
			# Compute the residuals
			e=Y-inprod.fdata(X.est,beta.est)
			
		}else if(type.basis=="pls"){
		
			if(is.null(p)) stop("Simple hypothesis with type.basis=\"pls\" need the number of components p")

			# Method
			meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in a PLS basis of ",p,"elements") 
			
			# Express X.fdata in a PLS basis
			fd2pls=fdata2plsr(X.fdata,Y,ncomp=p)
			if(length(fd2pls$l)!=1){
				X.est=fdata(mdata=fd2pls$x[,fd2pls$l]%*%fd2pls$rotation$data[fd2pls$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}else{
				X.est=fdata(mdata=fd2pls$x[,fd2pls$l]%*%t(fd2pls$rotation$data[fd2pls$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}		
			beta.est=beta0.fdata
			p.opt=p
			ord.opt=1:p.opt
			
			# Compute the residuals
			e=Y-inprod.fdata(X.est,beta.est)
			
		}else{
		
			stop(paste("Type of basis method",type.basis,"not implemented."))

		}
	
	}
	
	## 2. Bootstrap calibration ##
	
	# Start up
	pcvm.star=numeric(B)
	dif.beta.norms=numeric(B)
	
	# Calculus of the Aijr
	Aijr0=Aijr(X.est)
		
	# REAL WORLD
	pcvm=PCvM.statistic(X=X.est,residuals=e,p=p.opt,Aijr0=Aijr0)
	
	# BOOTSTRAP WORLD 
	if(show.prog) cat("Done.\nBootstrap calibration...\n ")
	if(show.prog) pb=txtProgressBar(style=3)
	for(i in 1:B){
		
		# Generate bootstrap errors
		e.hat=e*rber.gold(n)

		## COMPOSITE HYPOTHESIS ##
		if(is.null(beta0.fdata)){
		
			# Calculate Y.star
			Y.star=Y-e+e.hat
			
			# Calculate the estimated bootstrap errors
			if(est.method=="pc"){
			
				# Bootstrap model and residuals
				mod.pc.boot=fregre.pc.mod(fdataobj=X.fdata,y=Y.star,pc=pc.comp)
				e.hat.star[i,]=mod.pc.boot$residuals
				
				# Estimated bootstrap beta
				boot.beta.est[[i]]=mod.pc.boot$beta.est
				
			}else if(est.method=="pls"){
			
				# Bootstrap model and residuals
				mod.pls.boot=fregre.plsr.mod(fdataobj=X.fdata,y=Y.star,pls=pls.comp)
				e.hat.star[i,]=mod.pls.boot$residuals
				
				# Estimated bootstrap beta
				boot.beta.est[[i]]=mod.pls.boot$beta.est
			
			}else if(est.method=="basis"){
			
				# Bootstrap model and residuals
				mod.basis.boot=fregre.basis(fdataobj=X.fdata,y=Y.star,basis.b=basis.opt,basis.x=basis.opt)
				e.hat.star[i,]=mod.basis.boot$residuals
				
				# Estimated bootstrap beta
				boot.beta.est[[i]]=mod.basis.boot$beta.est
			
			}
			
			# Calculate PCVM.star
			pcvm.star[i]=PCvM.statistic(X=X.est,residuals=e.hat.star[i,],p=p.opt,Aijr0=Aijr0)
		
		## SIMPLE HYPOTHEIS ##
		}else{
		
			# Generate bootstrap errors
			e.hat.star[i,]=e.hat
		
			# Calculate PCVM.star
			pcvm.star[i]=PCvM.statistic(X=X.est,residuals=e.hat.star[i,],p=p.opt,Aijr0=Aijr0)
		
		}
		
		if(show.prog) setTxtProgressBar(pb,i/B)
		
	}

	## 3. MC estimation of the p-value and order the result ##
	
	# Compute the p-value
	pvalue<-sum(pcvm.star>pcvm)/B
	

	## 4. Graphical representation of the integrated process ##
	
	if(show.prog) cat("\nDone.\nComputing graphical representation... ")
	if(is.null(beta0.fdata) & plot.it){
	
		gamma=rproc2fdata(n=G,t=X.fdata$argvals)
		gamma=gamma/drop(norm.fdata(gamma))
		ind=drop(inprod.fdata(X.fdata,gamma))
		
		r=0.9*max(max(ind),-min(ind))
		u=seq(-r,r,l=200)
		mean.proc=numeric(length(u))
		mean.boot.proc=matrix(ncol=length(u),nrow=B.plot)
		
		res=apply(ind,2,function(iind){
			iind.sort=sort(iind,index.return=TRUE)
			stepfun(x=iind.sort$x,y=c(0,cumsum(e[iind.sort$ix])))(u)/sqrt(n)
		}
		)
		mean.proc=apply(res,1,mean)
	
		for(i in 1:B.plot){
		
			res=apply(ind,2,function(iind){
				iind.sort=sort(iind,index.return=TRUE)
				stepfun(x=iind.sort$x,y=c(0,cumsum(e.hat.star[i,iind.sort$ix])))(u)/sqrt(n)
			}
			)
			mean.boot.proc[i,]=apply(res,1,mean)
		
		}
		
		# Plot
		dev.new()
		main=ifelse(is.null(beta0.fdata),expression(paste("Integrated processs ",R[n]," for composite hypothesis")),expression(paste("Integrated processs ",R[n]," for simple hypothesis")))
		plot(u,mean.proc,ylim=c(min(mean.proc,mean.boot.proc),max(mean.proc,mean.boot.proc))*1.05,type="l",xlab=expression(paste(symbol("\341"),list(X, gamma),symbol("\361"))),ylab=expression(R[n](u)),main=main)
		for(i in 1:B.plot) lines(u,mean.boot.proc[i,],lty=2,col=gray(0.8))
		lines(u,mean.proc)
		text(x=0.75*u[1],y=0.75*min(mean.proc,mean.boot.proc),labels=sprintf("p-value=%.3f",pvalue))
		
	}
	if(show.prog) cat("Done.\n")

	
	# Result: class htest
	names(pcvm)="PCvM statistic"
	result<-structure(list(statistic=pcvm,boot.statistics=pcvm.star,p.value=pvalue,method=meth,B=B,type.basis=type.basis,
							beta.est=beta.est,boot.beta.est=boot.beta.est,p=p.opt,ord=ord.opt,data.name="Y=<X,b>+e"))
							
	class(result)<-"htest"
	return(result)
	
}
