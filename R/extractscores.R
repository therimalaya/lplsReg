#' Latent components extraction via NIPALS.
#'
#' Function for extracting latent components from L - structure data using the
#' NIPALS algorithm.
#'
#' @importFrom stats rnorm
#' @param X1 A data matrix (n x q) whose rows is connected to the rows of X2
#' @param X2 A data matrix (n x p) whose rows is connected to the rows of X1,
#' and whose columns are connected to the columns of X3.
#' @param X3 A data matrix (m x p) whose columns are connected toe the columns
#' of X2.
#' @param niter Number of NIPALS iterations.
#' @param threshold A threshold parameter in [0, 1) for potential
#' soft - thresholding of the rows of X3.
#' @return Returns a list og six latent vectors (normalized).
#' @author Solve Sæbø
#' @keywords latent lpls NIPALS
#' @examples
#'
#' data(BCdata)
#' scores <- extractscores(BCdata$Y,  BCdata$X,  BCdata$Z)
#' str(scores)
#' @export
extractscores <-
  function(X1, X2, X3, niter=5, threshold=0){
    t11 <- as.matrix(rnorm(dim(X1)[2], 0, 1))
    for (k in 1:niter) {
      t12 <- norm(projectonto(t(X1), t11))
      t22 <- norm(projectonto(X2, t12))
      if (is.null(X3)) {
        t31 <- t22
        t32 <- 0
      } else {
        t32 <- projectonto(t(X3), t22)
        if(threshold>0){
          t32 <- t32/max(abs(t32))
          t32 <- as.vector(sign(t32)*ifelse(abs(t32)<threshold, 0, abs(t32) - threshold))
        }
        t32 <- norm(t32)
        t31 <- norm(projectonto(X3, t32))
      }
      t21 <- norm(projectonto(t(X2), t31))
      t11 <- norm(projectonto(X1, t21))
    }
    scorelist <- list(t11=t11, t12=t12, t21=t21, t22=t22, t31=t31, t32=t32)
    return(scorelist)
  }
