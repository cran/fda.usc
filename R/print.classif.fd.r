print.classif.fd<-function (x, digits = max(3, getOption("digits") - 3), ...)
{
   cat("\n-Call:\n", deparse(x$C), "\n", sep = "")
#   cat("\n-Probability of correct classification by group (prob.classification):\n")
#   print(x$prob.classification)
if (x$C[[1]]=='classif.knn.fd'){
    cat("\n-Optimal number of neighbors: knn.opt=",x$knn.opt,
    "\nwith highest probability of correct classification max.prob=",
    x$max.prob,"\n")
    }
if (x$C[[1]]=='classif.kernel.fb'){
cat("\n-The greatest probability of correct classification: max.prob=",x$max.prob," is achieved with the number of bases: basis.opt=",x$basis.opt,"and  bandwidth value: h.opt=",x$h.opt,"\n")
 }
if (x$C[[1]]=='classif.kernel.fd'){
cat("\n-Optimal bandwidth: h.opt=",x$h.opt,"with highest probability of
correct classification: max.prob=",x$max.prob,"\n")
  }
    cat("\n")
    invisible(x)
}



