#################################################################
#################################################################
omit.fdata<-function(fdataobj,y=NULL){
nas<-apply(fdataobj$data,1,count.na)
if (!is.null(y)) {
nas.g<-is.na(y)
if (is.null(names(y))) names(y)<-1:length(y)
if (any(nas) & !any(nas.g)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
else {
if (!any(nas) & any(nas.g)) {
   cat("Warning: ",sum(nas.g)," values of group with NA are omited \n")
   bb<-!nas.g
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
else {
if (any(nas) & any(nas.g))  {
   bb<-!nas & !nas.g
   cat("Warning: ",sum(!bb)," curves  and values of group with NA are omited \n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
}}}
else {
if (any(nas)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
   }}
return(list(fdataobj,y))
}
#################################################################
#################################################################


missing.fdata<-function(fdataobj,basis=NULL){
   tt<-fdataobj$argvals
   rtt<-fdataobj$rangeval
   n<-nrow(fdataobj) 
   np<-length(tt)

   nas <- apply(fdataobj$data, 1, count.na)
   if (any(nas))  cat("Warning: ", sum(nas), " curves with NA are omited\n")
   nas <- which(nas)
   xall<-fdataobj
#   if (is.null(basis)) basis<-create.bspline.basis(rangeval = rtt, nbasis = max(5,floor(tt/4)))
   if (is.null(basis)) {
#      basis<-create.bspline.basis(rangeval = rtt, nbasis = length(tt))
      basis<-create.bspline.basis(rangeval = rtt, nbasis  = max(5,floor(tt/4)))
                       }
   for (i in nas) {
#  	cat(i)
	# is.na(xna$data[i,])
	curve<- which(!is.na(fdataobj$data[i,]))
	xall$data[i,-curve]<-eval.fd(tt,Data2fd(tt[curve],fdataobj$data[i,curve],basis))[-curve,1]
	}
xall
}

