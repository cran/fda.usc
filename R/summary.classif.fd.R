summary.classif.fd<-function (object, ...)
{
 cat("     - SUMMARY - \n")
 cat("\n-Probability of correct classification by group (prob.classification):\n")
 print(object$prob.classification)
 cat("\n-Confusion matrix between the theoretical groups (by rows)
  and estimated groups (by column) \n")
 print(table(object$group,object$group.est))
 if (object$C[[1]]=='classif.kernel.fb'){
  if (names(object[9])=="basis.opt") {
        cat("\n-Matrix of probability of correct classification
 by number of basis (num.basis) -rows- and bandwidth (h) -columns-:\n\n")
    print(round(1-object$misclassification,4))
        cat("\n-Functional measure of closeness (optimal distance, h.opt):\n")
         print(round(object$h.opt,4))
        cat("\n-Optimal number of basis: basis.opt=",object$basis.opt,"and optimal bandwidth: h.opt=",object$h.opt,
        "\nwith highest probability of correct classification: max.prob=",object$max.prob,"\n")
        
        cat("\n-The greatest probability of correct classification: max.prob=",object$max.prob,"
 is achieved with the number of bases: basis.opt=",object$basis.opt,"and  bandwidth value: h.opt=",object$h.opt,"\n")
    } }
 if (object$C[[1]]=='classif.kernel.fd'){
   cat("\n-Vector of probability of correct classification
    by banwidth (h):\n")
    print(round(1-object$misclassification,4))
   cat("\n-Functional measure of closeness (optimal distance, h.opt):\n")
   print(round(object$h.opt,4))

cat("\n-Optimal bandwidth: h.opt=",object$h.opt,"with highest probability of
correct classification: max.prob=",object$max.prob,"\n")
  }
 if (object$C[[1]]=='classif.knn.fd'){
   cat("\n-Vector of probability of correct classification
   by number of neighbors (knn):\n")
    print(round(1-object$misclassification,4))
    cat("\n-Optimal number of neighbors: knn.opt=",object$knn.opt,
    "\nwith highest probability of correct classification max.prob=",
    object$max.prob,"\n")
    }
cat("\n")
output<-object
}
