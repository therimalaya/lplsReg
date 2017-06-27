#' Error rate and confusion matrix
#'
#' The total error and the (g x g) confusion matrix for g-class classification
#' problems
#'
#' @importFrom stats model.matrix
#' @param trueclass The true classes
#' @param predclass The predicted classes
#' @param labs A vector of labels for the classes
#' @return
#' \item{total.error }{The total classification error}
#' \item{conf.matrix }{The (g x g) confusion matrix}
#' @author Solve Saebo
#' @keywords classification error
#' @examples
#'
#' data(BCdata)
#' confusion(factor(BCdata$Y), sample(c(0, 1), 130, replace = TRUE), c("0", "1"))
#' @export
confusion  <-  function(trueclass,  predclass,  labs){
  trueclass             <- factor(trueclass)
  predclass             <- factor(predclass)
  levels(predclass)     <- levels(trueclass)
  D1                    <- model.matrix(~trueclass - 1)
  D2                    <- model.matrix(~predclass - 1)
  conf.matrix           <- t(D1) %*% D2
  colnames(conf.matrix) <- rownames(conf.matrix) <-  labs
  temp                  <- conf.matrix
  diag(temp)            <- 0
  toterr                <- sum(temp)/length(trueclass)
  res                   <- (list(total.error = toterr, conf.matrix = conf.matrix))
  res
}
