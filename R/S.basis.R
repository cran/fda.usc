S.basis=function(tt,basis,lambda=0,Lfdobj=vec2Lfd(c(0,0)),...){
    phi=getbasismatrix(tt,basis)
    np<-length(tt)
    if (lambda!=0) {
       R=getbasispenalty(basis,Lfdobj)
       S=phi%*%solve(t(phi)%*%phi+lambda*R)%*%t(phi)}
   else {S=phi%*%solve(t(phi)%*%phi)%*%t(phi)}
    return(S)
}
