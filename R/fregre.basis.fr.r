fregre.basis.fr<- function(x,y,basis.s=NULL,basis.t=NULL,
lambda.s=0,lambda.t=0,Lfdobj.s=vec2Lfd(c(0,0),range.s),
Lfdobj.t=vec2Lfd(c(0,0),range.t),weights=NULL,...){
call<-match.call()
isfdx<-is.fd(x)
x.orig<-x
isfdy<-is.fd(y)
y.orig<-y
if (isfdx) {
  xfdobj<-x
  basis.x<- x$basis
  nbasis.x<- basis.x$nbasis
  range.x<- basis.x$rangeval
  if (is.null(basis.s))  basis.s<-basis.x
}
else {
  if (!is.fdata(x))  stop("x is not a functional data object of class fd or fdata.")
    range.x<-x$rangeval
    if (is.null(basis.s))  {
        np<-ncol(x)
        nbasis.b = min(51,floor(np/10))
        basis.s<-basis.x<-create.bspline.basis(rangeval=range.x,nbasis=nbasis.b)
   }
   else   basis.x<-basis.s
   xfdobj<-Data2fd(argvals =x$argvals, y = t(x$data), basisobj = basis.s,...)
}
if (isfdy) {
  yfdobj<-y
  basis.y<-y$basis
  nbasis.y<-basis.y$nbasis
  range.y<-basis.y$rangeval
#  if (!is.null(y$fdnames$time))  np.y<-length(y$fdnames$time)
#  else   np.y = max(c(51,2*nbasis.y+1))
  np.y<-nbasis.y
  tty <- seq(range.y[1],range.y[2],len=np.y)
  if (is.null(basis.t))  basis.t<-basis.y
}
else {
  if (!is.fdata(y))  stop("y is not a functional data object of class fd or fdata.")
  range.y<-y$rangeval
  tty<-y$argvals
  if (is.null(basis.t))  {
        np.y<-ncol(y)
        nbasis.y =max(min(51,floor(np.y/2)),7)
        basis.y<-basis.t<-create.bspline.basis(rangeval=range.y,nbasis=nbasis.y)
   }
   else   basis.y<-basis.t
   yfdobj<-Data2fd(argvals =y$argvals, y = t(y$data), basisobj = basis.t,...)
}

coefy   = yfdobj$coef
coefx   = xfdobj$coef

coefdx  = dim(coefx)
coefdy  = dim(coefy)
n=ncurves = coefdx[2]
################################################################################
alphabasis<-basis.t#create.constant.basis(range.y)
alphanbasis = alphabasis$nbasis

Finprod = inprod(basis.y, alphabasis)

alphattmat = diff(range.y) #cambiar
alphapenmat = NULL

################################################################################
range.t = basis.t$rangeval
range.s = basis.s$rangeval
if (range.s[1] != range.x[1] || range.s[2] != range.x[2]) {
    stop("Range of beta.s and x not compatible.")
}
nbasis.s = basis.s$nbasis
Hinprod = inprod(basis.x, basis.s)
xcoef = xfdobj$coef
basis.ss  = inprod(basis.s, basis.s)
S<-R<- NULL
################################################################################
range.t = basis.t$rangeval
if (range.t[1] != range.y[1] || range.t[2] != range.y[2]) {
    stop("Range of BETATFD coefficient and YFD not compatible.")
}
nbasis.t = basis.t$nbasis
Ginprod = inprod(basis.y, basis.t)
ycoef = yfdobj$coef
basis.tt  = inprod(basis.t, basis.t)


################################################################################
basis.talpha= inprod(basis.t, alphabasis)
################################################################################
Fmat = t(ycoef) %*% Finprod
Gmat = t(ycoef) %*% Ginprod
Hmat = t(xcoef) %*% Hinprod

if (is.null(weights)) {
    HHCP = t(Hmat) %*% Hmat
    HGCP = t(Hmat) %*% Gmat
    H1CP = as.matrix(colSums(Hmat))
    F1CP = as.matrix(colSums(Fmat))
} else {
    HHCP = t(Hmat) %*% (outer(weights,rep(nbasis.s))*Hmat)
    HGCP = t(Hmat) %*% (outer(weights,rep(nbasis.t))*Gmat)
    H1CP = t(Hmat) %*% weights
    F1CP = t(Fmat) %*% weights
}
################################################################################

betan = nbasis.s*nbasis.t
ncoef = alphanbasis+betan
Cmat  = matrix(0,ncoef,ncoef)
Dmat  = matrix(0,ncoef,1)

ind1 = 1:alphanbasis
ind2 = ind1
Cmat[ind1,ind2] = ncurves*alphattmat

ind2 = alphanbasis + (1:betan)
Cmat[ind1,ind2] = t(kronecker(H1CP,basis.talpha))
Dmat[ind1] = F1CP
ind1 =  alphanbasis  + (1:betan)
ind2 = 1:alphanbasis
Cmat[ind1,ind2] = t(Cmat[ind2,ind1])
ind2 = ind1
Cmat[ind1,ind2] = kronecker(HHCP,basis.tt)
if (lambda.s > 0) {
    R = eval.penalty(basis.s, Lfdobj.s)
    Cmat[ind1,ind2] = Cmat[ind1,ind2] + lambda.s*kronecker(R,basis.tt)
}
if (lambda.t > 0) {
    S = eval.penalty(basis.t, Lfdobj.t)
    Cmat[ind1,ind2] = Cmat[ind1,ind2] + lambda.t*kronecker(basis.ss,S)
}
Dmat[ind1] = matrix(t(HGCP),betan,1)
coefvec = qr.solve(Cmat, Dmat)
#coefvec = symsolve(Cmat, Dmat)
#  aa<-lm.fit(Cmat,Dmat)
basis.s$names
basis.t$names
ind1 = 1:alphanbasis
alpha.est<-coefvec[ind1]
beta.est<- matrix(coefvec[-ind1],nbasis.t,nbasis.s)
colnames(beta.est)<-basis.s$names
rownames(beta.est)<-basis.t$names
alpha.est= coefvec[ind1]
alphafdnames = yfdobj$fdnames
alphafdnames[[3]] = "Intercept"
alphafd = fd(alpha.est,  alphabasis, alphafdnames)
betafdnames = xfdobj$fdnames
betafdnames[[3]] = "Reg. Coefficient"                                        
betafd = bifd(t(beta.est), basis.s, basis.t, betafdnames)
beta.xest = beta.est %*% t(Hmat)
beta.xfd   = fd(beta.xest, basis.t)
yhat = eval.fd(tty, alphafd) %*% matrix(1,1,ncurves) + eval.fd(tty, beta.xfd)
if (isfdy) {
 fitted.values  <- fd(coef=yhat, basisobj=basis.y, fdnames=yfdobj$fdnames)
}
else {
 fitted.values<-fdata(t(yhat),y$argvals,y$rangeval,y$names)
}
residuals<-y-fitted.values
out = list(call=call,alpha.est=alphafd,coefficients=beta.est, beta.estbifd=betafd,
fitted.values=fitted.values,residuals=residuals,"lambda.s"=lambda.s,
"lambda.t"=lambda.t,Lfdobj.s=Lfdobj.s,Lfdobj.t=Lfdobj.t,weights=weights,
x=x,y=y,H=Hinprod,basis.s=basis.s,basis.t=basis.t,argvals.y=tty)
class(out)<-"fregre.fr"
return(out)
}
################################################################################


 #res03<-fregre.basis.fr(temp,prec,basis.t=smallbasis3,basis.s=smallbasis3) 