#' Function for lpls prediction.
#'
#' The function predicts new responses based on a fitted lpls-model and a set
#' of new predictor data.
#'
#'
#' @param object A fitted lpls-model as returned from \code{lplsReg}
#' @param ncomp The number of components to be used for prediction.
#' @param X2new A matrix of size (l x p) holding the predictor values of X2 for
#' l new observations.
#' @param silent Suppress warning about components being large
#' @return 
#' \item{pred}{The predicted response values.} 
#' \item{beta0 }{The estimated intercept in the lpls-regression model} 
#' \item{beta1}{The estimated vector of regression coefficients for the p 
#' predictor variables in X2.}
#' @author Solve S<c3><a6>b<c3><b8>
#' @keywords lpls prediction
#' @examples
#'
#' data(BCdata)
#' fit  <- lplsReg(BCdata$Y, BCdata$X, BCdata$Z, npc = 10)
#' pred <- predict(fit, 5, X2new = BCdata$X[1:10,])
#' pred$pred
#' @export
predict.lplsReg <- function(object, ncomp = object$npc,  X2new = NULL, silent = FALSE){

  reduced <- FALSE
  if (ncomp > object$maxcomp) {
      if (!silent) {
          cat("Warning: The selected component number is larger than in the fitted model. The largest component
          number in the fitted object will be used instead\n")
      }
    ncomp   <- object$maxcomp
    reduced <- TRUE
  }

  n     <- attr(object$var2, "dim")[1]
  p     <- attr(object$var2, "dim")[2]
  X2new <- as.matrix(X2new)
  ntest <- dim(X2new)[1]
  ptest <- dim(X2new)[2]

  if (sum(attr(object$var2, "rowm")) != 0) {
    rowm1 <- apply(X2new, 1, mean)
  } else {
    rowm1 <- rep(0, ntest)
  }

  beta1 <- object$T22[, 1:ncomp] %*%
    solve(t(object$P21[, 1:ncomp]) %*% object$T22[, 1:ncomp]) %*%
    t(object$P1[, 1:ncomp, drop = F])

  beta0 <- matrix(1, ntest, 1) %*%
    attr(object$var1, "colm") -
    (matrix(rep(1, ntest), ncol = 1) %*% attr(object$var2, "colm") +
     t(matrix(rep(1, p), ncol = 1) %*% rowm1) -
     matrix(attr(object$var2,"grandm"), nrow = ntest, ncol = p)) %*%
    beta1

  pred  <- beta0  +  X2new %*% beta1
  issig <- (beta1[, 1]!= 0) + 0
  res   <- list(pred = pred, beta0 = beta0, beta1 = beta1, reduced = reduced, issig = issig)
}
