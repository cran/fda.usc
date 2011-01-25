h.default=function (fdataobj, prob=c(0.025,0.25),len=51, metric = metric.lp,
                                        Ker = "AKer.norm", type.S ="S.NW",...)
{
    #type.S<-deparse(substitute(type.S))
    if (!is.fdata(fdataobj))   fdataobj = fdata(fdataobj)
    if (is.matrix(metric)) {mdist=metric}
    else {mdist=metric(fdataobj,fdataobj,...)}
    n=nrow(fdataobj)
#    x <- fdataobj[["data"]]
#    tt <- fdataobj[["argvals"]]
#    rtt <- fdataobj[["rangeval"]]
    if (type.S=="S.NW") {
#            mdist2 = mdist
            diag(mdist) = Inf
            h0 <- apply(mdist, 1, min, na.rm = TRUE)
            h.max = max(h0)
            h.med = median(h0)
            q.min = quantile(mdist, probs = prob[1], na.rm = TRUE)
            q.max = quantile(mdist, probs = prob[2], na.rm = TRUE)
            h.min = max(q.min, h.med)
            h.max = max(q.max, h.max)
            if (Ker== "AKer.norm") {
                    h.max = min(q.max, h.max)
                    h.min = min(q.min, h.med)
                  }
            h = unique(seq(h.min, h.max, len = len))
        }
        else if (type.S=="S.KNN") {
            if (len>n) {len=n-2;print("len=nrow-2")}
            h.min = floor(quantile(1:n, probs = prob[1], na.rm = TRUE,
                type = 4))
            h.max = floor(quantile(1:n, probs = prob[2], na.rm = TRUE,
                type = 4))
            h.min = max(2,h.min)
            h.max = max(h.min,h.max)
#            h =seq.int(h.min,h.max,len=len)
            h =unique(floor(seq(h.min,h.max,len=len)))
        }
        else {
        #           stop("Error in type.S argument
                }
        return(h)
 }

