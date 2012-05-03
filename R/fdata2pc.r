
################################################################
################################################################
fdata2pc<-function (fdataobj,  ncomp = 2,norm = TRUE,...)
{
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
    Jmin <- min(c(J, n))
    Xcen.fdata$data
    eigenres <- svd(Xcen.fdata$data)
    v <- eigenres$v
    u <- eigenres$u
    d <- eigenres$d
    D <- diag(d)
   delta=1
    vs <- fdata(t(v), tt, rtt, list(main = "fdata2pc", xlab = "t",
        ylab = "rotation"))
    scores <- matrix(0, ncol = J, nrow = n)
    if (norm) {
        dtt <- diff(tt)
        drtt <- diff(rtt)
        eps <- as.double(.Machine[[1]] * 10)
        inf <- dtt - eps
        sup <- dtt + eps
        if (all(dtt > inf) & all(dtt < sup))
            delta <- sqrt(drtt/(J - 1))
        else delta <- 1/sqrt(mean(1/dtt))
        no <- norm.fdata(vs)
        vs <- vs/delta
        newd <- d * delta
        scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs,...)
    }
    else {
        scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs, ...)
        newd <- d
    }
    colnames(scores) <- paste("PC", 1:J, sep = "")
    l <- 1:ncomp
    out <- list(call=C,d = newd, rotation = vs[1:ncomp],
         x = scores, fdataobj.cen = Xcen.fdata,
         mean = xmean, fdataobj = fdataobj,l=l)
    class(out) = "fdata.comp"
    return(out)
}

