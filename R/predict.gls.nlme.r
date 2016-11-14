predict.gls<-function (object, newdata, na.action = na.fail, ...)
{
#print("ppredict.gls.nlme.R predict.gls()")
    if (missing(newdata)) {
        return(fitted(object))
    }
#    print(object$model)
if (class(object$formula) == "formula")    form <- getCovariateFormula(object$formula)  
else   form <- getCovariateFormula(object)  
    mfArgs <- list(formula = form, data = newdata, na.action = na.action)
    mfArgs$drop.unused.levels <- TRUE
    dataMod <- do.call("model.frame", mfArgs)
    contr <- object$contrasts
    for (i in names(dataMod)) {
        if (inherits(dataMod[, i], "factor") && !is.null(contr[[i]])) {
            levs <- levels(dataMod[, i])
            levsC <- dimnames(contr[[i]])[[1]]
            if (any(wch <- is.na(match(levs, levsC)))) {
                stop(sprintf(ngettext(sum(wch), "level %s not allowed for %s",
                  "levels %s not allowed for %s"), paste(levs[wch],
                  collapse = ",")), domain = NA)
            }
            attr(dataMod[, i], "contrasts") <- contr[[i]][levs,
                , drop = FALSE]
        }
    }
    N <- nrow(dataMod)
    if (length(all.vars(form)) > 0) {
        X <- model.matrix(form, dataMod)
    }
    else {
        X <- array(1, c(N, 1), list(row.names(dataMod), "(Intercept)"))
    }
    cf <- coef(object)
    val <- c(X[, names(cf), drop = FALSE] %*% cf)
    attr(val, "label") <- "Predicted values"
    if (!is.null(aux <- attr(object, "units")$y)) {
        attr(val, "label") <- paste(attr(val, "label"), aux)
    }
    val
}
