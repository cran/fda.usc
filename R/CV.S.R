FCV.S <-function (y, S, W = diag(ncol(S)), trim = 0, draw = FALSE,metric=metric.lp,...)
{
    n = ncol(S)
    isfdata<-is.fdata(y)
    if (isfdata) {
       y2 = (y$data)     
       
       y.est = (S %*% y2)
#       digS<-diag(S)
#       I = diag(n)/(1 - digS)^2
#       W = W * I
       y.est<-fdata(y.est,y$argvals, y$rangeval, y$names)
       e <- y - y.est
#       ee <- drop(norm.fdata(e,metric=metric,...)[,1]^2)       
#       ee <- drop(norm.fdata(e,metric=metric,...)[,1]^2)
       ee<-sqrt(W)%*%(y2-y.est$data)
       ee<-fdata(ee,e$argvals,e$rangeval)            
       ee <- drop(norm.fdata(ee,metric=metric,...)[,1]^2)

       if (trim > 0) {
           e.trunc=quantile(ee,probs=(1-trim),na.rm=TRUE,type=4)
           ind <- ee <= e.trunc
           if (draw)     plot(y, col = (2 - ind))
           res = mean(ee[ind], na.rm = TRUE)
        }
        else  res = mean(ee, na.rm = TRUE)
    }
 if (is.nan(res))    res = Inf
 return(res)
}


CV.S <-function (y, S, W = diag(ncol(S)), trim = 0, draw = FALSE,metric=metric.lp,...)
{
    n = ncol(S)
    isfdata<-is.fdata(y)
    if (isfdata) { #representatation
       y2 = t(y$data)
       y.est = t(S %*% y2)
       y.est<-fdata(y.est,y$argvals, y$rangeval, y$names)
       e <- y - y.est
       ee <- norm.fdata(e,metric=metric,...)[,1]^2
       if (trim > 0) {
           e.trunc=quantile(ee,probs=(1-trim),na.rm=TRUE,type=4)
           ind <- ee <= e.trunc
           if (draw)     plot(y, col = (2 - ind))
           res = mean(ee[ind], na.rm = TRUE)
        }
        else  res = mean(ee, na.rm = TRUE)
    }
    else {
        y2 <- y
         y.est = S %*% y2
         I = diag(n)/(1 - diag(S))^2
         W = W * I
         e <- y2 - y.est
         if (trim > 0) {
            ee = t(e)
            e.trunc = quantile(abs(ee), probs = (1 - trim), na.rm = TRUE,
                type = 4)
            l <- which(abs(ee) <= e.trunc)
            res = mean(diag(W)[l] * e[l]^2, na.rm = TRUE)
         }
         res = mean(diag(W) * e^2)
    }
 if (is.nan(res))    res = Inf
 return(res)
}


