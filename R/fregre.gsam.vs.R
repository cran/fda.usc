################################################################################
fregre.gsam.vs <-function(data,y,x,alpha,
                      type.basis="pc",
                     #criterio="sp",
                      ncomp,
                      kbs,
                      #par.basis=list(),
                      dcor.min=.1,par.model,
                      xydist,
                      trace=FALSE,
                      CV=TRUE,
                      ncomp.fix=FALSE
                      ,smooth=TRUE
                      ){
  #falta insluir parámetro smooth!
  #0-
  criterio="sp"
  if (missing(y)) {stop("Hay que especificar el nombre de la variable respuesta")  }
  if (missing(alpha)) alpha<-.05
  if (missing(par.model)) par.model<-list() #se puede cambiar la familia
  n.edf.correction<-FALSE
  namdata<-names(data)
  idf<-which("df"==namdata )
  ydat<-data$df[,y]
  namfunc<-names(data[-idf])
  namnfunc<-setdiff(names(data$df),y)
#  print("***************     nombre variables *************")  
  #  print(namfunc)
  # print(namnfunc)
  #namnfunc<-setdiff(namnfunc,y)
   if (missing(x)) {  
     # print("misssing X")
     x<-c(namnfunc,namfunc)
     ifunc<-namfunc
     infunc<-namnfunc
   }
#     warning("Solo las variables funcionales son consideradas como covariables")
#   }
   else {
#print("hay Xss")    
    ifunc<-intersect(x,namfunc)
    infunc<-intersect(c(y,x),namnfunc)
#    xdatos<-as.list(data$df[,infunc,drop=F])
#    xdatos<-c(xdatos,data[ifunc])
   } 
# print("s1") 
 xdatos<-as.list(data$df[,y,drop=F])
 #print(names(xdatos)) 
#print("s2")  
 xdatos<-c(xdatos,as.list(data$df[,infunc,drop=F]),data[ifunc])
 
#print(names(xdatos)) 
#print("s3") 
# print(length(xdatos))    
  ldata0<-xdatos
  resp<-y
  xentran<-NULL
  tipoentran<-NULL
  #1- c?clulo distancias de cada objeto consigo mismo
#  print("1 calculando distancias")
  #    ldist0<-dist.list(ldata0)
  xynam<-names(ldata0)
  
  if (missing(xydist))    {
 #   print("s4")
    xydist<-ldist0<-dist.list(ldata0)
  #  print("s5")
  }    else {
    ldist0<-xydist[xynam]
  }
#  print("2 calculando distancias")
  #  print(names(ldist0))
  # print("fin names(ldist0")
  # dist_ij<-dcor.y(ldist0) 
  parar<-FALSE
  it<-1
  basisx<-NULL
  npredictors<-length(ldata0)
  ipredictors<-numeric(npredictors-1)
  names(ipredictors)<- setdiff(names(ldata0),resp)
  nvar<-npredictors-1
  
  if (missing(kbs)) kbs<- rep(-1,nvar)
  names(kbs)<-names(ipredictors)
  if (missing(ncomp)) ncomp<- rep(4,length(ifunc))
  names(ncomp)<-ifunc
  
  dcor=matrix(0,nrow= nvar,ncol= nvar)
  colnames(dcor)=c(names(ipredictors))
  rownames(dcor)=1:nvar
  #2- c?clulo correlaci?n  de cada la distancia de la respuesta vs distancia del resto de objetos
  # print("2 calculando correlaciones")
  # print(resp)
  # print(names(ldist0))
  dist_resp<-dcor.y(ldist0,resp)
  dcor[it,names(dist_resp)]<-dist_resp*(dist_resp > dcor.min)
  
  n.edf<-length(ydat)
  fpredictors.nl<-""
  form.nl<-paste(resp,"~",sep="")
  basis2<-list()  
  ycen<-data$df[,y]-mean(data$df[,y])
  gof<-NULL
  anyfdata<-FALSE
  while (!parar){
# print("3 Seleccion variable-Regresion")
    #      print("bucle")    
    esfactor<-FALSE #SE UTILIZA PARA NO PONER s(factor)
    #3- selecci?n dcor m?s elevada
    nam<-  names(dist_resp)
    ind.xentra<-which.max(dist_resp)
    xentra<-nam[ind.xentra]
    #dd<-dcor.ttest(ldist0[[resp]],ldist0[[xentra]],distance=TRUE)
    #print(n.edf)    
    #print(dd)
    dd<-dcor.test(ldist0[[resp]],ldist0[[xentra]],n=n.edf)
    if (trace){
      print(dd)
      print(xentra)
    }   
    if (par.fda.usc$verbose) {
      print(it)
      print(dist_resp)
      print(n.edf)
      print(nam)
      cat("Entra la variable ",xentra)
      print(dd)
      print(names(dd))
    }
    if (is.null(basisx)) {
      basisb<-basisx<-list()
    }
    if (dd$p.value>alpha) {
      parar=TRUE
      if (trace) print("the algorithm ends because no variable is significant")
    }    else{
      rownames(dcor)[it]<-xentra
      if (trace)      cat("Covariate: ")
      if (trace)       print(xentra)
      eslineal<-FALSE
      if (is.fdata(ldata0[[xentra]])) {
        # print(xentra)
        #if (xentra %in% names(par.basis)){
        #  itype<-which(xentra== names(par.basis))
        #  par.basis.ipred<-par.basis[[itype]]
        #}    else {par.basis.ipred <- list()}
        par.basis.ipred <- list()
        anyfdata<-TRUE
        par.basis.ipred$fdataobj<-ldata0[[xentra]]
        
        #if (is.list(type.basis))
        if (xentra %in% names(type.basis)){
          itype<-which(xentra== names(type.basis))
          type.basis.ipred<-type.basis[itype]
        }    else {type.basis.ipred <- type.basis[1]}
        
        
        tbasis<-type.basis.ipred
        if (type.basis.ipred=="fourier") tbasis<-"basis"
        if (type.basis.ipred=="bspline") tbasis<-"basis"
        
        nam<-"fregre.gsam.cv"
        #  par.basis.ipred<-list()
        par.basis.ipred$y<-resp#entra la etiqueta y no toda la variable como en el basis.cv o pc.cv
        par.basis.ipred$x<-xentra
        par.basis.ipred$data<-data

        # eslineal<-FALSE
        # if (linear) {
        #   res.linear<-flm.Ftest(data[[xentra]],data$df[[resp]], B = 1000, verbose = TRUE) 
        #   if (res.linear$p.value<.05) eslineal<-FALSE
        #   else eslineal <- TRUE
        # }

        # if (!eslineal)      
        #   res<- fregre.gsam.cv(data,resp,xentra, alpha=alpha,
        #                        type.basis=type.basis,criterio=criterio,
        #                        dcor.min=dcor.min,kmax=kmax)
        # if (eslineal){
        #      if (tbasis=="pc")
        #        res<-fregre.pc.cv(data[[xentra]],data$df[[resp]])
        #      if (tbasis=="basis")
        #        res<-fregre.basis.cv(data[[xentra]],data$df[[resp]])
        #      if (tbasis=="pls")
        #        res<-fregre.pls.cv(data[[xentra]],data$df[[resp]])             
        # }          
        if (CV)      {
          #print(type.basis)
          #print("CV")          
          res<- fregre.gsam.cv(data,resp,xentra, alpha=alpha,
                              type.basis=type.basis,
                              criterio=criterio,
                              ncomp=ncomp[xentra],ncomp.fix=ncomp.fix)  
          res$basis.x[[xentra]]$basis<-res$basis.x[[xentra]]$basis[res$pc.opt,]
#          res$basis.x[[xentra]]$x<-res$basis.x[[xentra]]$x[,res$pc.opt,drop=F]
          res$basis.x[[xentra]]$l<-res$pc.opt
          
          basisx[[xentra]]<-res$basis.x[[xentra]]
          basisb[[xentra]]<-res$basis.b[[xentra]]
         if (trace)   print("res fregre.gsam.cv")
         if (trace)   print(summary(res))
          
         # parar=TRUE
         }
        else{
          switch(tbasis,
         "pc"={
            best.pc<-1:ncomp[xentra]
            basis1<-create.pc.basis(ldata0[[xentra]],best.pc) #es lo mismo ?
          },
          "pls"={
            best.pc<-1:ncomp[xentra]
            basis1<-create.pls.basis(ldata0[[xentra]],ldata0[[resp]],best.pc)
          },"basis"={
#            best.pc<-1:kmax
            basis1<-create.bspline.basis(ldata0[[xentra]]$rangeval,nbasis=ncomp[xentra])
            basis2<-basis1
          })
          
          basisx[[xentra]]<-basis1
          basisb[[xentra]]<-basis2
        }
      }
      #     print(xentra)
      #     print(class(ldata0[[xentra]]))          
      if (!parar) {
        xentran<-c(xentran,xentra)
        ind.xentra2<-which(xentra==names(ldist0))
        ldist0<-ldist0[-ind.xentra2]
        ipredictors[xentra]<-ipredictors[xentra]+1
        #4- contrucci?n del modelo para esta variable #consido un cat?logo de 4 posibilidades (lineal/nolineal, funcional/scalar)
         if (is.factor(ldata0[[xentra]])) esfactor=TRUE                 
         #kbs2<-kbs
         fx=FALSE
           if (esfactor)      {
             fpredictors.nl<-xentra  
             # if (eslineal) { 
                   #   kbs=1
                   #  fx=TRUE
                   # }
             #print(fpredictors.nl)
                  }
        else         {    #hay factor el algun momento
            #print(" es factor o entra como lineal")
            #  print(fpredictors.nl)
            fpredictors.nl<-paste("s(",xentra,",k=",kbs,",","fx=",fx,")",sep="",collapse="+")
            
        }
#print("aa ver que pasa")        
#print(xentra);print(form.nl);print(fpredictors.nl)
         form.lin <-paste(form.nl,"+",xentra,sep="")
         form.nl<-(paste(form.nl,"+",fpredictors.nl,sep=""))  
        nam.model<-"fregre.gsam"          
        par.model$formula<-as.formula(form.nl)
        par.model$data<-data
        par.model[["basis.x"]]<-basisx
        res.nl<-do.call(nam.model,par.model)             
        if (!esfactor | !smooth){
          par.model.lin<-par.model
          par.model.lin$formula<-as.formula(form.lin )
          res.lin<-do.call(nam.model,par.model.lin)             
          a1=summary(res.nl)$sp.criterion
          a2=summary(res.lin)$sp.criterion
      #    print("iiiiiiiiiiiiiiiii")
          fact <-1.01
       #   cat(xentra,a1,a2,a2>(a1*fact),"\n")
          if (a2>(a1*fact)){
            res.nl<-res.lin
            form.1nl<-form.lin
            par.model<-par.model.lin
          }  	
        }
        if (trace)         print("sale fregre.gsam") #         }
       # if (trace) print(" ajuste del modelo GLS")
        #           print(nam.model)
        # print(names(par.model))          
        suma<-summary(res.nl)          
        dd <- res.nl$dims
        df <- dd[["p"]]
        edf<-dd[["N"]] -         dd[["p"]]    
        #         sr2 <- sum(res.nl$residuals^2)/edf
        r2 <- 1 - sum(res.nl$residuals^2)/sum(ycen^2)
        if (n.edf.correction) n.edf<- edf
        ldata0[[resp]]<-res.nl$residuals
        ldist0[[resp]]<-as.matrix(dist(ldata0[[resp]]), diag =TRUE, upper = TRUE,     p = 2)
        #2- c?clulo correlaci?n  de cada la distancia de la respuesta vs distancia del resto de objetos
#print(it)
#print(npredictors)
        if (it==(npredictors-1)) {
#print("iiiiittttttttt")          
          parar=TRUE
          #gof<-rbindgof,c(suma$logLik,suma$BIC,suma$AIC,edf,r2))
          gof<-rbind(gof,c(AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.sq,suma$dev.expl,res.nl$gcv.ubre))
        }
        else{
          it<-it+1
          #print("3 calculando correlaciones")
          dist_resp<-dcor.y(ldist0,resp)
          dcor[it,names(dist_resp)]=dist_resp*(dist_resp > dcor.min)
        }
        #     print(summary(res.nl))
        if (!parar){ 
          #            gof<-rbind(gof,drop(c(suma$logLik,suma$BIC,suma$AIC,edf,r2)))
          gof<-rbind(gof,c(AIC(res.nl),deviance(res.nl),
                           res.nl$df.residual,suma$r.sq,
                           suma$dev.expl,res.nl$gcv.ubre))
        }
      }  }
    
  }  
  # print(dim(gof));print(gof)
  gof<-data.frame(xentran,(gof))
  #  print(gof)
  #  print("aaa ver ohhh")
  #gof<-as.data.frame(gof)
  if (is.null(xentran)) {
    warning("ninguna variable seleccionada, se estima un LM a pelo")    
#print(as.formula(paste(resp,"~1",sep="")))    
#print(names(data))
    res.nl<-fregre.gsam(as.formula(paste(resp,"~1",sep="")),data=data)
    #res.nl<-gam(as.formula(paste(resp,"~1",sep="")),data=ldata0$df)
#print("peto")    
    suma<-summary(res.nl)
    gof<-data.frame(rbind(c(1,AIC(res.nl),deviance(res.nl),
                            res.nl$df.residual,suma$r.sq,
                            suma$dev.expl,res.nl$gcv.ubre)))      
  }  
  # print(gof)
  # print(dim(gof))"
  names(gof)<-c("xentra","AIC","deviance","df.residual","r.sq","dev.expl","GCV.ubre")
  out<-list("model"=res.nl,"gof"=gof,"i.predictor"=ipredictors,"xydist"=xydist
            ,"dcor"=dcor[1:(it-1),])
  class(out)<-c("select.gsam",class(out))
  return(out)
}
########################################################
