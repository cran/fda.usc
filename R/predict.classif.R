predict.classif<-function(object,new.fdataobj=NULL,type="class",...)
{
    lev<-levels(object$group)
    if (object$C[[1]] == "classif.glm") {
        prob <- ngroup <- length(table(object$group))
        if (ngroup == 2) {
            prob.group <- predict.fregre.glm(object$fit[[1]],
                newx = new.fdataobj, ...)
            yest <- ifelse(prob.group < 0.5,lev[1],lev[2])
        }
        else {
            prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]),
                ngroup))
            for (i in 1:ngroup) {
                prob.group[, i] <- predict.fregre.glm(object$fit[[i]],
                  newx = new.fdataobj, ...)
            }
            yest <- apply(prob.group, 1, which.min)
        }
        return(factor(yest,levels=lev))
    }
    if (!is.fdata(new.fdataobj))
        new.fdataobj = fdata(new.fdataobj, argvals(object$fdataobj))
    gg <- 1:nrow(new.fdataobj)
    nas <- apply(new.fdataobj$data, 1, count.na)
    if (any(nas)) {
        bb <- !nas
        cat("Warning: ", sum(nas), " curves with NA are omited\n")
        new.fdataobj$data <- new.fdataobj$data[bb, ]
        gg <- gg[bb]
    }
    tt = new.fdataobj[["argvals"]]
    new.fdata <- new.fdataobj[["data"]]
    mdatos = object$fdataobj[["data"]]
    nn = ncol(mdatos)
    numgr = nrow(mdatos)
    numgrn = nrow(new.fdata)
    if (object$C[[1]] == "classif.kernel") {
        as <- list()
        C <- object$C
        m <- object$m
        imetric <- m[5]
        if (imetric == 0) {
            a1 <- metric.lp
            len.metricc <- length(formals(a1))
            vv <- array(0, dim = c(len.metricc))
        }
        else {
            a1 <- match.fun(C[[imetric]])
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
                }
                else {
                  as[[ind.m]] <- formals(a1)[[ind.m]]
                }
                ind.m <- ind.m + 1
            }
        }
        if (is.vector(new.fdata)) {
            print("1 new data are predicted")
            if (length(new.fdata) != dim(mdatos)[2]) {
                stop("ERROR IN THE DATA DIMENSIONS")
            }
            nndatos = array(NA, dim = c(1, nn))
            numgrn = 1
            nndatos[1, ] = new.fdata
            tr = 1
        }
        else {
            if (is.matrix(new.fdata)) {
                if (dim(mdatos)[2] != dim(new.fdata)[2]) {
                  stop("ERROR IN THE DATA DIMENSIONS")
                }
                else {
                  numgrn = dim(new.fdata)[1]
                  cat(numgrn, " new data are predicted  \n")
                  nndatos = new.fdata
                  tr = 3
                }
            }
            else {
                stop("ERROR IN THE DATA DIMENSIONS")
            }
        }
        if ((dim(mdatos)[1]) != (length(object$group))) {
            stop("ERROR IN THE DATA DIMENSIONS")
        }
        group = object$group
        if (!is.factor(group)) {
            group = as.factor(group)
        }
        numg = nlevels(group)
        lev<-levels(group)
        group.pred = array(0, dim = numgrn)
        prob.groupV2 = array(0, dim = c(numg))
        pr = 1
        maxk = 0
        prob = 0
        mm = mdatos
        mdist = array(0, dim = c(numgr, numgrn))
        pr = 1
        prob = 0
        vec = 0
        as[[2]] <- new.fdataobj
        as[[1]] <- object$fdataobj
        mdist <- do.call(a1, as)
        vmax = array(0, dim = c(numgr))
        vvmax = array(0, dim = c(numgr))
        mat = array(0, dim = c(numgr, numg))
        hmdist = mdist/object$h.opt
        pgrup = array(0, dim = c(numg, numgrn))
        for (j in 1:numg) {
            y = as.integer(group == levels(group)[j])
            pgrup[j, ] <- apply(hmdist, 2, rkernel, y = y)
        }
        l = array(0, dim = c(numgrn))
        for (i in 1:numgrn) {
            l[i] = which.max(pgrup[, i])
            group.pred[i] = levels(group)[l[i]]
        }
    group.pred=factor(group.pred,levels=lev)
    if (type=="class")   return(group.pred)
    else {
        colnames(pgrup) <- rownames(new.fdataobj$data)
        rownames(pgrup) <- levels(object$group)
        return(list(group.pred = group.pred, prob.group = pgrup))
    }
    }
    else {
        if (object$C[[1]] == "classif.knn") {
            as <- list()
            C <- object$C
            m <- object$m
            imetric <- m[5]
            if (imetric == 0) {
                a1 <- metric.lp
                len.metricc <- length(formals(a1))
                vv <- array(0, dim = c(len.metricc))
            }
            else {
                a1 <- match.fun(C[[imetric]])
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
                  }
                  else {
                    as[[ind.m]] <- formals(a1)[[ind.m]]
                  }
                  ind.m <- ind.m + 1
                }
            }
            if (is.vector(new.fdata)) {
                print("1 new data are predicted")
                if (length(new.fdata) != dim(mdatos)[2]) {
                  stop("ERROR IN THE DATA DIMENSIONS")
                }
                nndatos = array(NA, dim = c(1, nn))
                numgrn = 1
                nndatos[1, ] = new.fdata
                tr = 1
            }
            else {
                if (is.matrix(new.fdata)) {
                  if (dim(mdatos)[2] != dim(new.fdata)[2]) {
                    stop("ERROR IN THE DATA DIMENSIONS")
                  }
                  else {
                    numgrn = dim(new.fdata)[1]
                    cat(numgrn, "new data are predicted \n")
                    nndatos = array(NA, dim = c(dim(new.fdata)[1],
                      dim(new.fdata)[2]))
                    nndatos = new.fdata
                    tr = 3
                  }
                }
                else {
                  stop("DATA ERROR")
                }
            }
            group<-object$group
            if (!is.factor(object$group)) group = as.factor(group)
            numg = nlevels(group)
            group.pred = array(0, dim = c(numgrn))
            dist = array(0, dim = c(numgr, numgrn))
            dist2 = array(0, dim = c(numgr, numgrn))
            pr = 1
            prob = 0
            vec = array(0, dim = c(numgrn))
            as[[2]] <- new.fdataobj
            as[[1]] <- object$fdataobj
            dist <- do.call(a1, as)
            pgrup = array(0, dim = c(numg, numgrn))
            for (k in 1:numgrn) {
                vvmax = vmax = array(0, dim = c(numgr))
                vec = quantile(dist[, k], prob = (object$knn.opt/(numgr)),
                  type = 4)
                l = which(dist[, k] <= vec)
                l2 = which.min(dist[, k])
                for (j in 1:object$knn.opt) pgrup[group[l[j]],
                  k] = pgrup[group[l[j]], k] + 1
                max2 = which.max(pgrup[, k])
                if (length(max2) > 1) {
                  vvmax = as.character(group[l2])
                }
                else {
                  vvmax = levels(group)[max2]
                }
                group.pred[k] = vvmax
            }
            pgrup <- pgrup/object$knn.opt
        }
    }
    group.pred=factor(group.pred,levels=lev)
    if (type=="class")   return(group.pred)
    else {
        colnames(pgrup) <- rownames(new.fdataobj$data)
        rownames(pgrup) <- levels(object$group)
        return(list(group.pred = group.pred, prob.group = pgrup))
    }
}



