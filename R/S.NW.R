S.NW<-function (tt, h, Ker = Ker.norm,cv=FALSE)
{
 if (is.matrix(tt)) {
    if (ncol(tt)!=nrow(tt)) {
      if (ncol(tt)==1) {
         tt=as.vector(tt)
         tt=abs(outer(tt,tt, "-"))}
      #else stop("Error: incorrect arguments passed")
    }}
 else if (is.vector(tt))    tt=abs(outer(tt,tt, "-"))
 else stop("Error: incorrect arguments passed")
  if (cv)  diag(tt)=Inf
  k=Ker(tt/h)
#  if (cv)   diag(k)=0
  S =k/apply(k,1,sum)

return(S)}


#res.np=fregre.np.cv(x,y,Ker=AKer.tri,type.CV=CV.S)
#res.np$gcv

# porque sino habría que calcular el cuantil de I*e^2 que es el efecto conjunto del error y el peso

#S.NW pasar cv
#truncamiento?
