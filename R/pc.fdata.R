pc.fdata=function(out,l=NULL,loadings=FALSE,draw=TRUE,...) {
  if (class(out)=="fregre.fd") {
    if (!is.null(l)) l=l
    else l=out$l
   dat<-out$fdataobj$data
   pr.com=out$svd.fdata
   out=pr.com$x
   rotation=aperm(pr.com$rotation$data)
 }
 else {
  if (is.null(l)) l=1:3
  if (!is.fdata(out)) stop("No fdata class")
 #   pr.com2=prcomp(out,center=TRUE) ###
   dat=out$data
   pr.com=pc.svd.fdata(out) ###
   rotation=aperm(pr.com$rotation$data)
    }
 le=length(l)
 pc=rep(NA,le)
 if (loadings) {  ##
    pr.com=varimax(rotation,normalize = FALSE)
    p=pr.com$loadings
 }
else p=pr.com$x
pr.x= apply(p, 2, var)/sum(var(p))
 if (draw){
  C<-match.call()
  lenC=length(C)
  j=1
  while (j<=lenC) {
        if (names(C)[j]=="ask") {
           ask=C[[j]]
           j=lenC +1
           }
        else { j=j+1; ask=FALSE  }
  }
   dev.new()
   if (ask) {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          for (i in 1:le) {
          ts.plot(rotation[,l[i]],ylab=c("loadings",l[i],sep=""),
      main=c(paste("PC",l[i],"- Expl. Var. ",round(pr.x[l[i]] * 100, 2),"%",sep="")))
      if (i<le)
      for (j in (i+1):le) {

            if (nrow(out)<50)   {
                        plot(p[,c(i,l[j])],main="BIPLOT",type="n")
                        text(p[,c(i,l[j])],rownames(out))
                        }
             else                         plot(p[,c(i,l[j])],main="BIPLOT")
            if (nrow(out)<50)      {
                           plot(p[,c(l[j],i)],main="BIPLOT",type="n")
                           text(p[,c(l[j],i)],rownames(out))
               }
           else  plot(p[,c(l[j],i)],main="BIPLOT")
        } }  }
    else   {
    par(mfrow=c(le,le))
    for (i in 1:le) {
      par(mfg=c(i,i))
      ts.plot(rotation[,l[i]],ylab=c("loadings",l[i],sep=""),
       main=c(paste("PC",l[i],"- Expl. Var. ",round(pr.x[l[i]] * 100, 2),"%",sep="")))
      if (i<le)
      for (j in (i+1):le) {
            par(mfg=c(i,j))
            if (nrow(out)<50)     {
                      plot(p[,c(i,l[j])],main="BIPLOT")
                      text(p[,c(i,l[j])],rownames(out))
                      }
            else plot(p[,c(i,l[j])],main="BIPLOT")

            par(mfg=c(j,i))
            if (nrow(out)<50)            {
               plot(p[,c(l[j],i)],main="BIPLOT",type="n")
               text(p[,c(l[j],i)],rownames(out))
               }
            else plot(p[,c(l[j],i)],main="BIPLOT")
     }  }    }  }
  cat("\n-With",le,"Principal Components is  explained ",round(sum(pr.x[l])*100
 ,2),"%\n of the variability of explicative variables. \n
-Variability for each  principal components -PC- (%):\n")
print(round(pr.x[l] * 100, 2))
 return(invisible(pr.com))
}

