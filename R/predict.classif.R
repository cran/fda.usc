predict.classif<-function(object,new.fdataobj=NULL,type="class",...)
{
    if (is.null(new.fdataobj)) return(object$group.est)
    isfdata<-is.fdata(new.fdataobj)
    object$group<-factor(object$group,levels=levels(object$group)[which(table(object$group)>0)])
    ny<-lev<-levels(object$group)
#    if (is.null(object)) stop("No classif object entered")
#    if (is.null(new.fdataobj)) stop("No newx entered")
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
                 group.pred<-(factor(yest,levels=lev))
                 pgrup<-prob.group #######
    }
  if (object$C[[1]] == "classif.gsam") {
        prob <- ngroup <- length(table(object$group))
        if (ngroup == 2) {
            prob.group <- predict.fregre.gsam(object$fit[[1]],
                newx = new.fdataobj, ...)
            yest <- ifelse(prob.group < 0.5,lev[1],lev[2])
        }
        else {
            prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]),
                ngroup))
            for (i in 1:ngroup) {
                prob.group[, i] <- predict.fregre.gsam(object$fit[[i]],
                  newx = new.fdataobj, ...)
            }
            yest <- apply(prob.group, 1, which.min)
        }
                 group.pred<-(factor(yest,levels=lev))
                  pgrup<-prob.group #######
    }
    if (object$C[[1]] == "classif.gkam") {
        prob <- ngroup <- length(table(object$group))
        if (ngroup == 2) {

            prob.group <- predict.fregre.gkam(object$fit[[1]],
               newx = new.fdataobj,...)
            yest <- ifelse(prob.group < 0.5,lev[1],lev[2])
        }
        else {
            prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]),
                ngroup))
            for (i in 1:ngroup) {
                prob.group[, i] <- predict.fregre.gkam(object$fit[[i]],
                  newx = new.fdataobj, ...)
            }
            yest <- apply(prob.group, 1, which.min)
        }
                 group.pred<-(factor(yest,levels=lev))
                  pgrup<-prob.group #######
    }
   
if (object$C[[1]] == "classif.np") {
#if (!is.fdata(new.fdataobj)) new.fdataobj=fdata(new.fdataobj,object$fdataobj[["argvals"]],object$fdataobj[["rangeval"]],object$fdataobj[["names"]])
gg<-1:nrow(new.fdataobj)
if (isfdata) {
nas<-apply(new.fdataobj$data,1,count.na)
if (any(nas)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   new.fdataobj$data<-new.fdataobj$data[bb,]
   gg<-gg[bb]
   }
  newx<-new.fdataobj[["data"]]
  tt<-new.fdataobj[["argvals"]]
  rtt<-new.fdataobj[["rangeval"]]
}
else newx<-as.matrix(new.fdataobj)
nn <- nrow(new.fdataobj)
if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
# if (is.vector(newx))  newx <- t(as.matrix(newx))
 x=object$fdataobj
 y=object$y
  h=object$h.opt
# h<-0.5
 n = nrow(x)
 nn = nrow(newx)
 np <- ncol(x)
 ny<-levels(y)
 numg=nlevels(y)
# if (n != (length(y)))         stop("ERROR IN THE DATA DIMENSIONS")
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
 bs=as<-list()
 C<-object$call
 m<-object$m
 Ker=object$Ker
   par.metric<-list()
   par.metric<-attr(object$mdist,"par.metric")
   parm<-attr(object$mdist,"par.metric")
   a1<-attr(object$mdist,"call")
if (isfdata) {
   par.metric[["fdata2"]]<-x
   par.metric[["fdata1"]]<-new.fdataobj
   nmdist <- do.call(a1,par.metric)
   }
else {
   par.metric[["x"]]<-new.fdataobj
   par.metric[["y"]]<-x
   nmdist <- (do.call(a1,par.metric))
  }
#  print(nmdist)
   object$par.S$cv=FALSE
   object$par.S$tt<-nmdist
   kmdist=object$type.S(nmdist,h=h,Ker=object$Ker, cv = FALSE)
         pgrup = array(0, dim = c(numg, nn))
        l = array(0, dim = c(nn))
        group.pred = array(0, dim = nn)
        for (j in 1:numg) {
            grup = as.integer(y == lev[j])
            pgrup[j, ] <- kmdist%*%matrix(grup,ncol=1)
        }
    group.pred<-factor(ny[apply(pgrup,2,which.max)],levels=ny)
#    group.pred<-factor(ifelse(apply(pgrup,2,which.max)==1,lev[1],lev[2])   )
#return(group.pred)
}
#    group.pred=factor(group.pred,levels=lev)

if (type=="class")   return(group.pred)
else {
    if (isfdata)          colnames(pgrup) <- rownames(new.fdataobj$data)
    else           colnames(pgrup) <- rownames(new.fdataobj)
        rownames(pgrup) <- levels(object$group)
        return(list(group.pred = group.pred, prob.group = pgrup))
    }
}


