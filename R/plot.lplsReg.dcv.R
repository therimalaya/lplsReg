#' Double-CV plot
#'
#' Plot function summarizing some results from double-crossvalidation of
#' lpls-regression/classification.
#'
#' The first plot is a plot of posterior probability of class membership for
#' each sample plotted versus sample number. For g-group classification, there
#' will be g dots per sample, and the largest dot indicates the predicted
#' class/group. The second plot is a plot of the frequency of each variable
#' being selected by jack-knifing in each of the cross-validation segments.
#' High selection frequency may be considered as a measure of variable
#' importance and stability.
#'
#' @param object A double-Cv object as returned from \code{lplsReg.dcv}
#' @param identifyVariable Logical. Should interactive variable identification
#' be activated?
#' @author Solve Saebo
#' @keywords posterior variable-selection
#' @export
plot.lplsReg.dcv <- function(object, identifyVariable = FALSE){
  dev.new()
  matplot(
    x    = 1:length(object$apost[, 1]),
    y    = object$apost,
    pch  = 20,
    cex  = 1.5,
    col  = object$trueclass,
    main = "Posterior probabilities for class membership",
    xlab = "Sample number",
    ylab = "Posterior probability"
  )
  dev.new()
  plot(
    x    = 1:length(object$varfreq),
    y    = object$varfreq,
    type = "h",
    main = "Frequency of variables selected in Jackknifing",
    xlab = "Variable",
    ylab = "Frequency"
  )
  selected <- NULL
  if (identifyVariable) {
    cat("Use the mouse to identify variables\n")
    selected <- identify(1:length(object$varfreq), object$varfreq)
  }
  print(selected)
}
