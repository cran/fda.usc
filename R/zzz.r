globalVariables('icnt')



#' @keywords internal

# .onLoad <- function(libname, pkgname){
#     #.libPaths("~/RLibrary")
#     require(foreach)
#     require(parallel)
#     require(doParallel)
#     ncores = ops.fda.usc()$ncores
#     cat("--------------------------------------------------------------\n")
#     cat("Functional Data Analysis in R \n")
#    #  cat("Loaded fda.usc 2.0-1\n")
#     cat("fda.usc version 2.0-0 (built on 2019-11-11) is now loaded\n")
#     cat("fda.usc will use",ncores,"cores in foreach function\n")
#     cat("ops.fda.usc function changes the parameters of the package\n")
#     #cat("Please, use ops.fda.usc function to change global package parameters\n")
#     cat("--------------------------------------------------------------\n")
# 	NULL
#   }
#   
#   
# .onAttach <- function (lib, pkg) {
#   pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="geoR"),
#                             fields=c("Title","Version","Date")))
#   
#   packageStartupMessage(paste("--------------------------------------------------------------\n",
#                               pkg.info["Title"]),"\n",
#                         " For an Introduction to geoR go to http://www.leg.ufpr.br/geoR\n",
#                         paste(" geoR version ", pkg.info["Version"],
#                               " (built on ", pkg.info["Date"], ") is now loaded\n", sep=""),
#                         "--------------------------------------------------------------\n"
#   )
# }

  # .onLoad(pkgname="fda.usc")
  
  #.onAttach <- function (lib, pkg) {
#   .onLoad <- function(libname, pkgname){ 
#     #require(foreach);    require(parallel);    require(doParallel)
#     
#     requireNamespace("foreach", quietly = TRUE)
#     requireNamespace("parallel", quietly = TRUE)
#     requireNamespace("doParallel", quietly = TRUE)
#     #loadNamespace("parallel");    print(1);    loadNamespace("doParallel");    print(2)  ;    loadNamespace("foreach")    print(3);    #loadNamespace("fda.usc")
#       ncores = max(parallel::detectCores() -1,1)-3
#     print(4)    
#     .par.fda.usc = list()
#     .par.fda.usc$verbose = FALSE
#     .par.fda.usc$trace = FALSE
#     .par.fda.usc$warning = FALSE
#     .par.fda.usc$ncores = ncores
#     .par.fda.usc$int.method = "TRAPZ"
#     .par.fda.usc$eps =  as.double(.Machine[[1]]*10)
#     print(44)  
#     if (ncores==1) {
#       registerDoSEQ()
#       print(55)
#     }  else{
#       print(66)
#       print(ncores)
#         #cl <-  parallel::makePSOCKcluster(ncores )
#         print(77)
#         #doParallel::registerDoParallel(cl)
#     }
#     print(88)
#     e<-environment(ops.fda.usc)
#     print(e)
#     unlockBinding("par.fda.usc", e)
#     assign("par.fda.usc", .par.fda.usc, envir = e)
#     get("par.fda.usc", envir = e)
#     
#     print(99)
#     pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="fda.usc"),
#                               fields=c("Title","Version","Date")))
#     print(101)
#    packageStartupMessage(
#      paste("-----------------------------------------------------------------\n",pkg.info["Title"]),"\n",
# 		#	 "Functional Data Analysis in R \n",
#     paste(" fda.usc version ", pkg.info["Version"]," (built on ", pkg.info["Date"], ") is now loaded\n", sep=""),
# 			  paste(" fda.usc will use",ncores,"cores in foreach function\n"),
# 			  " ops.fda.usc function changes the parameters of the package\n",
# 			 "-----------------------------------------------------------------\n"
#     )
#    print(111)
#    return(invisible())
# }

  
#' @import fda 
#' @import splines
#' @import MASS
#' @import mgcv
#' @import stats
#' @import nlme
  
  # @import rpart deleted
  
  #import(foreach,"getDoParWorkers","getDoParRegistered","getDoParName","getDoParVersion","foreach","registerDoSEQ","setDoSeq","setDoPar")
  #importFrom(parallel, "makeCluster", "stopCluster", "detectCores", "clusterExport", "clusterEvalQ")
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators  icount

#' @importFrom foreach %dopar% getDoParWorkers getDoParRegistered getDoParName getDoParVersion foreach registerDoSEQ setDoSeq setDoPar
#' @import parallel
#' @importFrom grDevices adjustcolor colorRampPalette dev.cur dev.interactive dev.new devAskNewPage gray heat.colors palette
#' @importFrom graphics abline boxplot contour curve filled.contour image legend lines matlines pairs par persp plot points polygon rect rug stars text title
#' @importFrom methods callGeneric
#' @importFrom utils combn installed.packages modifyList setTxtProgressBar txtProgressBar globalVariables
  
  
 # .onLoad <- function(libname, pkgname){ 
 #    #require(foreach);    require(parallel);    require(doParallel)
 #    
 #    requireNamespace("foreach", quietly = TRUE)
 #    requireNamespace("parallel", quietly = TRUE)
 #    requireNamespace("doParallel", quietly = TRUE)
 #    #loadNamespace("parallel");    print(1);    loadNamespace("doParallel");    print(2)  ;    loadNamespace("foreach")    print(3);    #loadNamespace("fda.usc")
 #    # setHook(packageEvent("fda.usc","attach"),function(pkgname, libpath) {
 #    #   ncores <- max(parallel::detectCores() -1,1)
 #    #    if (ncores==1) {
 #    #      foreach::registerDoSEQ()
 #    #    }  else{
 #    #      if (foreach::getDoParWorkers()!=ncores){
 #    #        # cat("getDoParWorkers != ncores")
 #    #        cl <-  suppressWarnings(parallel::makePSOCKcluster(ncores ))
 #    #        doParallel::registerDoParallel(cl)
 #    #      }
 #    #    }
 #    # }    )
 #    
 #    ncores <- max(parallel::detectCores() -1,1)
 #    .par.fda.usc<-list()
 #    .par.fda.usc$verbose = FALSE
 #    .par.fda.usc$trace = FALSE
 #    .par.fda.usc$warning = FALSE
 #    .par.fda.usc$ncores = ncores
 #    .par.fda.usc$int.method = "TRAPZ"
 #    .par.fda.usc$eps = as.double(.Machine[[1]]*10)
 #    e<-environment(ops.fda.usc)
 #    unlockBinding("par.fda.usc", e)
 #    assign("par.fda.usc", .par.fda.usc, envir = e)
 #    get("par.fda.usc", envir = e)
 #    #if (is.null(ncores)) 
 #    #   ncores = max(parallel::detectCores() -1,1)
 #    # #print("entra ops.fda.usc")
 #    # .par.fda.usc = list()
 #    # .par.fda.usc$verbose = FALSE
 #    # .par.fda.usc$trace = FALSE
 #    # .par.fda.usc$warning = FALSE
 #    # .par.fda.usc$ncores = ncores
 #    # .par.fda.usc$int.method = "TRAPZ"
 #    # .par.fda.usc$eps =  as.double(.Machine[[1]]*10)
 #    # if (ncores==1) {
 #    #   foreach::registerDoSEQ()
 #    # }  else{
 #    #   if (foreach::getDoParWorkers()!=ncores){
 #    #     # cat("getDoParWorkers != ncores")
 #    #     cl <-  suppressWarnings(parallel::makePSOCKcluster(ncores ))
 #    #     doParallel::registerDoParallel(cl)
 #    #   }
 #    # }
 #    return(invisible())
 #  }

  .onAttach <- function(lib, pkg,...){
    pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="fda.usc"),
                              fields=c("Title","Version","Date")))
    foo <- suppressWarnings(foreach::"%dopar%"(foreach::foreach(i=1), {}))
     packageStartupMessage(
      paste("-----------------------------------------------------------------\n",pkg.info["Title"]),"\n",
      #	 "Functional Data Analysis in R \n",
      paste(" fda.usc version ", pkg.info["Version"]," (built on ", pkg.info["Date"], ") is now loaded\n", sep=""),
      paste(" fda.usc is running sequentially usign foreach package\n"),
      paste(" Please, execute ops.fda.usc() once to run in local parallel mode\n"),
      #" ops.fda.usc() changes the parameters of the package\n",
      "-----------------------------------------------------------------\n"
    )
  }
  
  # in gbmCrossVal.R
  #  gbmCluster<- function (n) 
  # {
  #   if (is.null(n)) {
  #     n <- parallel::detectCores()
  #   }
  #   parallel::makeCluster(n)
  # }
  #Set up cluster and add finalizer
  # cluster <- gbmCluster(n.cores)
  # on.exit(parallel::stopCluster(cluster))

  
  #fdasrvf:::kmeans.R
  #if (parallel){
  #  cores <- detectCores()
  #  cl <- makeCluster(cores)
  # registerDoParallel(cl)
    #} else
    #{
    #  registerDoSEQ()
  #}  
  #if (parallel){
   # stopCluster(cl)
  #}