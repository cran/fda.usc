
##############################################
##     Auxiliary functions for flm.test   	##
##############################################

##############################################
## File created by Eduardo García-Portugués ##
## using code from library fda.usc          ##
##############################################

# Modification of the fda.usc fregre.basis.cv function to force the same number of elements in the basis of x and beta
fregre.basis.cv.mod=function(x,y,p=seq(5,45,by=1),type.basis="bspline",type.CV=GCV.S,silent=TRUE,...){

	res=numeric()
	for(i in 1:length(p)){
		if(!silent) cat(p[i],"\n")
		res[i]=fregre.basis.cv(x,y,basis.b=p[i],basis.x=p[i],type.basis=type.basis,type.CV=type.CV,...)$gcv.opt
	}
	opt=p[which.min(res)]
	basis=eval(parse(text=paste("create.",type.basis,".basis(rangeval=x$rangeval,nbasis=opt)",sep="")))

	return(fregre.basis(x,y,basis.b=basis,basis.x=basis))

}

# Modification of the fregre.pc function to use the previously computed PC without recomputing
fregre.pc.mod=function(fdataobj,y,pc){
    if (!is.fdata(fdataobj)) 
        fdataobj = fdata(fdataobj)
    nas <- apply(fdataobj$data, 1, count.na)
    nas.g <- is.na(y)
    if (is.null(names(y))) 
        names(y) <- 1:length(y)
    if (any(nas) & !any(nas.g)) {
        bb <- !nas
        cat("Warning: ", sum(nas), " curves with NA are omited\n")
        fdataobj$data <- fdataobj$data[bb, ]
        y <- y[bb]
    }
    else {
        if (!any(nas) & any(nas.g)) {
            cat("Warning: ", sum(nas.g), " values of group with NA are omited \n")
            bb <- !nas.g
            fdataobj$data <- fdataobj$data[bb, ]
            y <- y[bb]
        }
        else {
            if (any(nas) & any(nas.g)) {
                bb <- !nas & !nas.g
                cat("Warning: ", sum(!bb), " curves  and values of group with NA are omited \n")
                fdataobj$data <- fdataobj$data[bb, ]
                y <- y[bb]
            }
        }
    }
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x)
    np <- ncol(x)
	l <- pc$l # NEW
    lenl = length(l)
    if (n != (length(y))) 
        stop("ERROR IN THE DATA DIMENSIONS")
    C <- match.call()
    if (is.null(rownames(x))) 
        rownames(x) <- 1:n
    ycen = y - mean(y)
	# NEW pc <- fdata2pc(fdataobj, norm, ncomp = max(l))
    xcen <- pc$fdataobj.cen
    if (length(l) == 1) {
        vs <- pc$rotation$data[l, ]
        Z <- pc$x[, l]
    }
    else {
        vs <- t(pc$rotation$data[l, ])
        Z <- (pc$x[, l])
    }
    cnames <- colnames(pc$x)[l]
    response = "y"
    df <- data.frame(y, Z)
    colnames(df) <- c("y", cnames)
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf, "+", cnames[i], 
        sep = "")
    object.lm = lm(formula = pf, data = df, x = TRUE, y = TRUE)
    beta.est <- object.lm$coefficients[2:(lenl + 1)] * pc$rotation[l, 
        ]
    beta.est$data <- apply(beta.est$data, 2, sum)
    beta.est$names$main <- "beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data), nrow = 1)
    Z = cbind(rep(1, len = n), Z)
    S = solve(t(Z) %*% Z)
    H <- Z %*% S %*% t(Z)
    e <- object.lm$residuals
    df = lenl + 1
    sr2 <- sum(e^2)/(n - df)
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    out <- list(call = C, beta.est = beta.est, fitted.values = object.lm$fitted.values, 
        fdata.comp = pc, coefficients = object.lm$coefficients, 
        residuals = object.lm$residuals, df = df, r2 = r2, sr2 = sr2, 
        H = H, fdataobj = fdataobj, y = y, l = l, lm = object.lm, 
        pc = pc)
    class(out) = "fregre.fd"
    return(out)
}

# Modification of the fregre.pls function to use the previously computed PLS without recomputing
fregre.plsr.mod=function(fdataobj,y,pls){
    if (!is.fdata(fdataobj)) 
        fdataobj = fdata(fdataobj)
    nas <- apply(fdataobj$data, 1, count.na)
    nas.g <- is.na(y)
    if (is.null(names(y))) 
        names(y) <- 1:length(y)
    if (any(nas) & !any(nas.g)) {
        bb <- !nas
        cat("Warning: ", sum(nas), " curves with NA are omited\n")
        fdataobj$data <- fdataobj$data[bb, ]
        y <- y[bb]
    }
    else {
        if (!any(nas) & any(nas.g)) {
            cat("Warning: ", sum(nas.g), " values of group with NA are omited \n")
            bb <- !nas.g
            fdataobj$data <- fdataobj$data[bb, ]
            y <- y[bb]
        }
        else {
            if (any(nas) & any(nas.g)) {
                bb <- !nas & !nas.g
                cat("Warning: ", sum(!bb), " curves  and values of group with NA are omited \n")
                fdataobj$data <- fdataobj$data[bb, ]
                y <- y[bb]
            }
        }
    }
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x)
    np <- ncol(x)
	pc <- pls # NEW
	l <- pls$l # NEW
    lenl = length(l)
    if (n != (length(y))) 
        stop("ERROR IN THE DATA DIMENSIONS")
    C <- match.call()
    if (is.null(rownames(x))) 
        rownames(x) <- 1:n
    ycen = y - mean(y)
    # NEW pc <- fdata2plsr(fdataobj, ycen, ncomp = max(l))
    xcen <- pc$fdataobj.cen
    if (length(l) == 1) {
        vs <- pc$rotation$data[l, ]
        Z <- pc$x[, 1:l]
    }
    else {
        vs <- pc$rotation$data[l, ]
        Z <- (pc$x[, l])
    }
    cnames <- colnames(pc$x)[l]
    response = "y"
    df <- data.frame(y, Z)
    colnames(df) <- c("y", cnames)
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf, "+", cnames[i], 
        sep = "")
    object.lm = lm(formula = pf, data = df, x = TRUE, y = TRUE)
    beta.est <- object.lm$coefficients[2:(lenl + 1)] * pc$rotation[l, 
        ]
    beta.est$data <- apply(beta.est$data, 2, sum)
    beta.est$names$main <- "beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data), nrow = 1)
    Z = cbind(rep(1, len = n), Z)
    S = solve(t(Z) %*% Z)
    H <- Z %*% S %*% t(Z)
    e <- object.lm$residuals
    df = lenl + 1
    sr2 <- sum(e^2)/(n - df)
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    out <- list(call = C, beta.est = beta.est, fitted.values = object.lm$fitted.values, 
        fdata.comp = pc, coefficients = object.lm$coefficients, 
        residuals = object.lm$residuals, df = df, r2 = r2, sr2 = sr2, 
        H = H, fdataobj = fdataobj, y = y, l = l, lm = object.lm)
    class(out) = "fregre.fd"
    return(out)
}

# Modification of the min.basis function in order to not show screen information
min.basis.mod=function(fdataobj,type.CV=GCV.S,W=diag(1,ncol=np,nrow=np),lambda=0,numbasis=floor(seq(ncol(fdataobj)/16,ncol(fdataobj)/2,len=10)),type.basis="bspline",par.CV=list(trim=0,draw=FALSE),...) {
    if (!is.fdata(fdataobj)) 
        fdataobj = fdata(fdataobj)
    nas1 <- apply(fdataobj$data, 1, count.na)
    if (any(nas1)) 
        stop("fdataobj contain ", sum(nas1), " curves with some NA value \n")
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    nam <- fdataobj[["nam"]]
    lenlambda <- length(lambda)
    lenbasis <- length(numbasis)
    nc <- nrow(fdataobj)
    np <- ncol(fdataobj)
    gcv <- array(Inf, dim = c(lenbasis, lenlambda))
    GCV.basis.min = Inf
    as <- list()
    as[[1]] <- rtt
    names(as)[[1]] <- "rangeval"
    C <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("fdataobj", "tt", "type.CV", "W", "lambda", 
        "numbasis", "type.basis"), names(mf), 0L)
    imetric <- m[7]
    if (imetric == 0) {
        a1 <- create.bspline.basis
        len.metricc <- length(formals(a1))
        vv <- array(0, dim = c(len.metricc))
    }
    else {
        a1 <- paste("create.", type.basis, ".basis", sep = "")
        len.metricc <- length(formals(a1))
        vv <- array(0, dim = c(len.metricc))
    }
    ii <- imetric + 1
    if (C[ii] != "NULL()") {
        ind.m <- 3
        while (C[ii] != "NULL()" && ind.m <= len.metricc) {
            aa <- any(names(C) == names(formals(a1))[ind.m])
            if (aa) {
                vv[ind.m] <- which(names(C) == names(formals(a1)[ind.m]))
                ii <- ii + 1
                as[[ind.m]] <- C[[vv[ind.m]]]
                names(as)[[ind.m]] <- names(formals(a1)[ind.m])
            }
            else {
                as[[ind.m]] <- formals(a1)[[ind.m]]
            }
            ind.m <- ind.m + 1
        }
    }
    for (i in 1:lenbasis) {
        as[[2]] <- numbasis[i]
        names(as)[[2]] <- "nbasis"
        base <- do.call(a1, as)
        for (k in 1:lenlambda) {
            S2 <- S.basis(tt, base, lambda[k])
            gcv[i, k] <- type.CV(fdataobj, S = S2, W = W, trim = par.CV$trim, 
                draw = par.CV$draw, ...)
        }
    }
    l = which.min(gcv)
    i = (l%%lenbasis)
    k = (l%/%lenbasis) + 1
    if (i == 0) {
        i = lenbasis
        k = k - 1
    }
    lambda.opt <- lambda[k]
    numbasis.opt <- numbasis[i]
    gcv.opt <- gcv[l]
    as[[2]] = numbasis[i]
    names(as)[[2]] <- "nbasis"
    base.opt = do.call(a1, as)
    S.opt <- S.basis(tt, base.opt, lambda[k])
    fdata.est <- S.opt %*% t(x)
    # cat("\n The minimum GCV (GCV.OPT=", round(gcv.opt, 4), sep = "", 
        # ") is achieved with\n the number of basis (numbasis.opt=", 
        # numbasis.opt, ")\n and lambda value    (lambda.opt=", 
        # lambda.opt, ")\n\n")
    if (length(numbasis) > 1) 
        dimnames(gcv)[[1]] <- numbasis
    if (length(lambda) > 1) 
        dimnames(gcv)[[2]] <- lambda
    if (lenbasis > 1) {
        if (numbasis.opt == min(numbasis)) 
            cat(" Warning: numbasis.opt is the minimum number of basis provided, range(numbasis)=", 
                range(numbasis), "\n")
        else if (numbasis.opt == max(numbasis)) 
            cat(" Warning: numbasis.opt is the maximum number of basis provided, range(numbasis)=", 
                range(numbasis), "\n")
    }
    if (lenlambda > 1) {
        if (lambda.opt == min(lambda)) 
            cat(" Warning: lambda.opt is the minimum lambda value provided, range(lambda)=", 
                range(lambda), "\n")
        else if (lambda.opt == max(lambda)) 
            cat(" Warning: lambda.opt is the maximum lambda value provided, range(lambda)=", 
                range(lambda), "\n")
    }
    fdata.est = fdata(t(fdata.est), tt, rtt, nam)
    output <- list(gcv = gcv, numbasis = numbasis, lambda = lambda, 
        fdataobj = fdataobj, fdata.est = fdata.est, gcv.opt = gcv.opt, 
        numbasis.opt = numbasis.opt, lambda.opt = lambda.opt, 
        S.opt = S.opt, base.opt = base.opt)
}

# Old fregre.pls of fda.usc that uses the method of Preda (2002)
fregre.plsr=function(fdataobj,y,l=1:3,...){
    if (!is.fdata(fdataobj)) 
        fdataobj = fdata(fdataobj)
    nas <- apply(fdataobj$data, 1, count.na)
    nas.g <- is.na(y)
    if (is.null(names(y))) 
        names(y) <- 1:length(y)
    if (any(nas) & !any(nas.g)) {
        bb <- !nas
        cat("Warning: ", sum(nas), " curves with NA are omited\n")
        fdataobj$data <- fdataobj$data[bb, ]
        y <- y[bb]
    }
    else {
        if (!any(nas) & any(nas.g)) {
            cat("Warning: ", sum(nas.g), " values of group with NA are omited \n")
            bb <- !nas.g
            fdataobj$data <- fdataobj$data[bb, ]
            y <- y[bb]
        }
        else {
            if (any(nas) & any(nas.g)) {
                bb <- !nas & !nas.g
                cat("Warning: ", sum(!bb), " curves  and values of group with NA are omited \n")
                fdataobj$data <- fdataobj$data[bb, ]
                y <- y[bb]
            }
        }
    }
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    names <- fdataobj[["names"]]
    n = nrow(x)
    np <- ncol(x)
    lenl = length(l)
    if (n != (length(y))) 
        stop("ERROR IN THE DATA DIMENSIONS")
    C <- match.call()
    if (is.null(rownames(x))) 
        rownames(x) <- 1:n
    ycen = y - mean(y)
    pc <- fdata2plsr(fdataobj, ycen, ncomp = max(l))
    xcen <- pc$fdataobj.cen
    if (length(l) == 1) {
        vs <- pc$rotation$data[l, ]
        Z <- pc$x[, 1:l]
    }
    else {
        vs <- pc$rotation$data[l, ]
        Z <- (pc$x[, l])
    }
    cnames <- colnames(pc$x)[l]
    response = "y"
    df <- data.frame(y, Z)
    colnames(df) <- c("y", cnames)
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf, "+", cnames[i], 
        sep = "")
    object.lm = lm(formula = pf, data = df, x = TRUE, y = TRUE)
    beta.est <- object.lm$coefficients[2:(lenl + 1)] * pc$rotation[l, 
        ]
    beta.est$data <- apply(beta.est$data, 2, sum)
    beta.est$names$main <- "beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data), nrow = 1)
    Z = cbind(rep(1, len = n), Z)
    S = solve(t(Z) %*% Z)
    H <- Z %*% S %*% t(Z)
    e <- object.lm$residuals
    df = lenl + 1
    sr2 <- sum(e^2)/(n - df)
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    out <- list(call = C, beta.est = beta.est, fitted.values = object.lm$fitted.values, 
        fdata.comp = pc, coefficients = object.lm$coefficients, 
        residuals = object.lm$residuals, df = df, r2 = r2, sr2 = sr2, 
        H = H, fdataobj = fdataobj, y = y, l = l, lm = object.lm)
    class(out) = "fregre.fd"
    return(out)
}

# Old fregre.pls.cv of fda.usc that uses the method of Preda (2002)
fregre.plsr.cv=function (fdataobj, y, kmax = 8, criteria = "CV", ...){
    if (!is.fdata(fdataobj)) 
        fdataobj = fdata(fdataobj)
    nas <- apply(fdataobj$data, 1, count.na)
    nas.g <- is.na(y)
    if (is.null(names(y))) 
        names(y) <- 1:length(y)
    if (any(nas) & !any(nas.g)) {
        bb <- !nas
        cat("Warning: ", sum(nas), " curves with NA are omited\n")
        fdataobj$data <- fdataobj$data[bb, ]
        y <- y[bb]
    }
    else {
        if (!any(nas) & any(nas.g)) {
            cat("Warning: ", sum(nas.g), " values of group with NA are omited \n")
            bb <- !nas.g
            fdataobj$data <- fdataobj$data[bb, ]
            y <- y[bb]
        }
        else {
            if (any(nas) & any(nas.g)) {
                bb <- !nas & !nas.g
                cat("Warning: ", sum(!bb), " curves  and values of group with NA are omited \n")
                fdataobj$data <- fdataobj$data[bb, ]
                y <- y[bb]
            }
        }
    }
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    n <- nrow(x)
    nc <- ncol(x)
    ind = 1:kmax
    l = l2 = list()
    ck = 1
    tab = list("AIC", "AICc", "SIC", "SICc", "HQIC", "rho", "CV")
    type.i = pmatch(criteria, tab)
    if (is.na(type.i)) 
        stop("Error: incorrect criteria")
    else {
        if (type.i < 7) {
            cv.AIC <- rep(NA, kmax)
            for (j in 1:kmax) {
                out = fregre.plsr(fdataobj, y, l = 1:ind[j], ...)
                ck <- length(l) + 1
                s2 <- sum(out$residuals^2)/n
                ck = ind[j] + 1
                if (criteria == "AIC") {
                  cv.AIC[j] <- log(s2) + 2 * (ck)/n
                }
                else if (criteria == "AICc") {
                  cv.AIC[j] <- log(s2) + 2 * (ck)/(n - ck - 2)
                }
                else if (criteria == "SIC") {
                  cv.AIC[j] <- log(s2) + log(n) * ck/n
                }
                else if (criteria == "SICc") {
                  cv.AIC[j] <- log(s2) + log(n) * ck/(n - ck - 
                    2)
                }
                else if (criteria == "HQIC") {
                  cv.AIC[j] <- log(s2) + 2 * log(log(n)) * ck/n
                }
                else if (criteria == "rho") {
                  A <- out$residuals
                  B <- 1 - diag(out$H)/n
                  D1 <- (A/B)^2
                  cv.AIC[j] <- sum(D1)
                }
            }
            min.AIC = min(cv.AIC)
            pc.opt <- 1:ind[which.min(cv.AIC)]
        }
        else {
            cv.AIC <- rep(NA, kmax)
            out = fregre.plsr(fdataobj, y, l = 1:kmax, ...)
            cv.AIC <- out$fdata.comp$res.pls$validation$PRESS
            min.AIC = min(cv.AIC)
            pc.opt <- 1:which.min(cv.AIC)
        }
    }
    names(cv.AIC) = paste("PLS", 1:kmax, sep = "")
    fregre = fregre.plsr(fdataobj, y, l = pc.opt, ...)
    return(list(fregre.pls = fregre, pls.opt = pc.opt, MSC.min = min.AIC, 
        MSC = cv.AIC))
}

# Old fdata2pls of fda.usc that uses the method of Preda (2002)
fdata2plsr=function (fdataobj, y, ncomp = 2, norm = TRUE, ...){
    C <- match.call()
    if (!is.fdata(fdataobj)) 
        stop("No fdata class")
    nas1 <- apply(fdataobj$data, 1, count.na)
    if (any(nas1)) 
        stop("fdataobj contain ", sum(nas1), " curves with some NA value \n")
    X <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
    nam <- fdataobj[["names"]]
    mm <- fdata.cen(fdataobj)
    xmean <- mm$meanX
    Xcen.fdata <- mm$Xcen
    n <- nrow(Xcen.fdata)
    J <- ncol(Xcen.fdata)
    Jmin <- min(c(ncomp, J, n))
    cnames <- paste("PLS", seq(1, ncol(fdataobj$data)), sep = "")
    ycen = y - mean(y)
    response = "y"
    df <- data.frame(y, Xcen.fdata$data)
    colnames(df) <- c("y", cnames)
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) pf <- paste(pf, "+", cnames[i], 
        sep = "")
    res <- plsr(as.formula(pf), ncomp = ncomp, data = df, method = "oscorespls", 
        validation = "LO", ...)
    class(res$loadings) <- "matrix"
    vs <- fdata(t(res$loadings), tt, rtt, list(main = "pls.fdata", 
        xlab = "t", ylab = "rotationloadings"))
    scores <- matrix(0, ncol = J, nrow = n)
    if (norm) {
        no <- norm.fdata(vs)
        vs$data <- sweep(vs$data, 1, drop(no), "/")
        scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs)
    }
    else {
        scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs)
    }
    colnames(scores) <- paste("PLS", 1:J, sep = "")
    out <- list(rotation = vs, x = scores, res.pls = res, fdataobj.cen = Xcen.fdata, 
        mean = xmean, y = y, l = 1:ncomp, C = C)
    class(out) = "fdata.comp"
    return(out)
}

