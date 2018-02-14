predict.fregre.igls<-function (object, new.fdataobj = NULL, data, 
          se.fit = FALSE, scale = NULL, df = df, interval = "none",
          level = 0.95, weights = 1, pred.var, n.ahead =1L,...) 
{
  if (is.null(object)) 
    stop("No object entered")
  if (is.null(new.fdataobj)) {
    yp <- object$fitted.values
    print("No newx entered")
    return(yp)
  }
  else {
    yp <- predict.fregre.lm(object, new.fdataobj, se.fit = se.fit, 
                            scale = scale, df = df, interval = interval, level = level, 
                            weights = weights, pred.var = pred.var, ...)
    predictor <- yp
    nn <- length(yp)
    n <- length(object$residuals)
    predictor <- drop(yp)
    res.var <- object$sr2
    if (missing(pred.var)) pred.var = res.var/weights
    if (!is.null(object$corStruct)) {
      if (names(object$correlation) == "cor.AR" | names(object$correlation) == 
          "cor.ARMA") {
          if ((class(object$corStruct[[1]])[1] == "Arima" | 
             class(object$corStruct[[1]])[1] == "ar") & 
            length(object$corStruct) > 1) {
          ype <- NULL
          gr <- data[, object[["correlation"]][[1]][["group"]]]
          lev <- levels(gr)
          tab <- table(gr)
          for (j in 1:length(tab)) {
            lennn <- tab[j]
            previousone <- object$corStruct[[lev[j]]]
            ind <- gr == lev[j]
            if (lennn != 0) {
              if (class(object$corStruct[[1]])[1] == 
                  "Arima") 
                ype[ind] = predict(object$corStruct[[lev[j]]], 
                                   se.fit = se.fit, n.ahead = lennn)
            }
          }
          }
        if (class(object$corStruct$ar) == "Arima") {
          ype = predict(object$corStruct$ar, se.fit = se.fit, 
                        n.ahead = nn)
          }
        if (class(object$corStruct$ar) == "ar") {
          ype = predict(object$corStruct$ar, object$residuals, 
                        se.fit = se.fit, n.ahead = nn)
          }          
        if (class(object$corStruct[[1]])[1] == "lm") {
          coef.lm <- coef(object$corStruct$lm)
          p <- length(coef.lm)
          lenp <- length(p)
          ype <- NULL
          gr <- data[, object[["correlation"]][[1]][["group"]]]
          if (!is.factor(gr)) 
            gr <- factor(gr)
          lev <- levels(gr)
          tab <- table(gr)
          for (j in 1:length(tab)) {
            lennn <- tab[j]
            previousone <- object$corStruct$lm$res.x[, 
                                                     j]
            ind <- gr == lev[j]
            e <- rep(NA, len = lennn)
            for (i in 1:lennn) {
              e[i] <- sum(coef.lm * previousone)
              previousone <- c(e[i], previousone[-p])
            }
            ype[ind] <- e
          }
        }
      }
      predictor <- predictor + as.numeric(ype)
    }
    return(predictor)
  }
  return(predictor)
}