#' Summary function for lplsReg double cross-validation objects.
#'
#' Summary function for lplsReg double cross-validation objects.
#'
#'
#' @param obj An object returned from a call to \code{lplsReg.dcv}.
#' @return Prints the call, the total classification error and the confusion
#' matrix.
#' @author Solve S<c3><a6>b<c3><b8>
#' @keywords classification error-rate
#' @export
summary.lplsReg.dcv <- function(obj){
  cat("\nCall:\n", paste(deparse(obj$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("\nTotal classification error:\n",obj$total.error[[1]], sep = " ")
  cat("\n\nConfusion matrix:\n")
  print(obj$total.error[[2]])
}
