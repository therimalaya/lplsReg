#' LPLS - model regression
#'
#' Based on the NIPALS algorithm this function fits an LPLS - regression model to
#' a dummy - response vector X1 using a predictor matrix X2 with potentially
#' background knowledge matrix X3. This function only supports a single
#' response variable,  that is,  X1 must be a vector. The fitted model is similar
#' to a regular PLS - model in case X3 is NULL. The matrix X3 typically hold
#' background information about the variables in X2 (e.g. root matrix of some
#' prior covariance matrix,  or relevant old data). Based on the size of the
#' parameter alpha,  the background information may alter the regular
#' PLS - loading weights in order to make the background information influential
#' on the estimated regression coefficients of the final linear model.
#'
#' @importFrom stats model.matrix cor
#' @param X1 A response vector or matrix for regression. For classification
#' this should be either a factor or a dummy coded 0/1 matrix with one column
#' per group.
#' @param X2 Predictor matrix of size (n x p).
#' @param X3 Background information matrix of size (m x p)
#' @param npc Number of components to use in the LPLS model.
#' @param alpha Parameter between 0 and 1 for controlling the influence of X3.
#' If 0 no influence,  if 1 maximum influence.
#' @param rowcentering A logical vector of length 2. The elements setting the
#' potential row - centering of X2 and X3,  respectively.  Default is c(FALSE,
#' FALSE),  that is,  no row centering.
#' @param colcentering A logical vector of length 2. The elements setting the
#' potential column centering of X2 and X3,  respectively.  Default is c(TRUE,
#' FALSE),  that is,  column centering of X2 only. If X3 holds a set of \code{m}
#' extra observations of the same type as X2,  but with no obsereved response,
#' then X3 and X2 should be centered in the same manner.
#' @param grandcentering A logical vector of length 2. If TRUE the respective
#' matrix (X2 or X3) will be centered using the overall mean only.
#' @param pathshrink A numeric [0, 1). May be used to perform soft - shrinkage of
#' the rows of X3 as a type of variable selection on the background
#' information. The properties of this is so far not much explored. See S<c3><a6>b<c3><b8>
#' et al. 2008b on ST - PLS for a similar approach for variable selection in PLS.
#' @param niter The number of NIPALS iterations. Default is 10.
#' @return
#' \item{call}{The function call used}
#' \item{npc}{The number of components used}
#' \item{P1}{Correlation loadings vector for X1}
#' \item{P21}{Correlation loadings for X2}
#' \item{P22}{Correlation scores for X2}
#' \item{P3}{Correlation loadings for X3}
#' \item{T11 - T32 }{Latent vector matrices for
#' X1,  X2 and X3 (a mix of loading weights and scores in PLS terminology).}
#' \item{P1 }{Loading vectors for the variables in X1}
#' \item{P21 }{Loading vectors for the variables in X2}
#' \item{P3 }{Loading vectors for X3}
#' \item{var1 - var3 }{The (centered) data matrices used as input}
#' \item{var1res - var3res }{Residual matrices after extracting the information
#' from the \code{npc} latent components.}
#' @author Solve S<c3><a6>b<c3><b8>
#' @references S<c3><a6>b<c3><b8>,  S.,  Alm<c3><b8>y,  T.,  Flatberg,  A.,  Aastveit,  A.H.,  Martens,  H.
#' (2008a) LPLS - regression: a method for prediction and classification under
#' the influence of background information on predictor variables. Chemometrics
#' and Intelligent Laboratory Systems,  91 (2) 121 - 132.
#'
#' S<c3><a6>b<c3><b8>,  S.,  Alm<c3><b8>y,  T.,  Aaroe,  J.,  Aastveit,  A.H. (2008b) ST - PLS: A
#' multi - directional nearest shrunken centroid type classifier via Partial
#' Least Squares. Journal of Chemometrics,  22 (1),  54 - 62.
#' @keywords lpls regression
#' @examples
#'
#' data(BCdata)
#' fit.class <- lplsReg(factor(BCdata$Y),  BCdata$X,  BCdata$Z,  npc = 10)
#' plot(fit.class)
#' #For regression,  drop the factor() statement and treat the response as a continuous variable.
#' fit.reg <- lplsReg(BCdata$Y,  BCdata$X,  BCdata$Z,  npc = 10)
#' @export
lplsReg <-
  function(X1,  X2,  X3,  npc = 2,  alpha = 0.5,
           rowcentering = c(F, F),  colcentering = c(T, F), grandcentering = c(F, F),
           pathshrink = 0,  niter = 30){

    #Dimensions
    if(class(X1) == "factor"){
      X1 <- model.matrix(~X1 - 1)
    }
    X1dim <- dim(X1)
    X2dim <- dim(X2)
    if(!is.null(X3)){
        X3dim <- dim(X3)
        }else{
        X3dim <- c(1, X2dim[2])
    }

    #Centering
    X1 <- centering(X1, FALSE, TRUE, FALSE)
    X2 <- centering(X2, rowcentering[1], colcentering[1], grandcentering[1])
    if(!is.null(X3)){
        X3 <- centering(X3, rowcentering[2], colcentering[2], grandcentering[2])
    }

#-------------------------------------------------------------

    locX1 <- X1
    locX2 <- X2
    locX3 <- X3

    T11 <- matrix(0, nrow = X1dim[2], ncol = npc)
    T12 <- matrix(0, nrow = X1dim[1], ncol = npc)
    T21 <- matrix(0, nrow = X2dim[1], ncol = npc)
    T22 <- matrix(0, nrow = X2dim[2], ncol = npc)
    T31 <- matrix(0, nrow = X3dim[2], ncol = npc)
    T32 <- matrix(0, nrow = X3dim[1], ncol = npc)
    P1  <- matrix(0, nrow = X1dim[2], ncol = npc)
    P3  <- matrix(0, nrow = X3dim[1], ncol = npc)
    P21 <- matrix(0, nrow = X2dim[2], ncol = npc)
    for(a in 1:npc){
        latent <- extractscores(locX1, locX2, locX3, niter = niter, threshold = pathshrink)

        #Eigenvectors and scores
        T11[, a] <- latent$t11
        T12[, a] <- latent$t12
        T21[, a] <- latent$t21
        T22[, a] <- latent$t22
        T31[, a] <- latent$t31
        T32[, a] <- latent$t32

        #A new weight - vector t22 is constructed as a intermediate version of original t22 and t31
        T22[, a] <- norm(alpha*T31[, a] + (1 - alpha) * T22[, a])
        T21[, a] <- locX2 %*% T22[, a]

        #Estimation of regression vectors
        P1[, a]                  <- t(locX1) %*% T21[, a] %*% solve(crossprod(T21[, a]))
        if(!is.null(X3)){P3[, a] <- locX3 %*% T22[, a] %*% solve(crossprod(T22[, a]))}
        P21[, a]                 <- t(locX2) %*% T21[, a] %*% solve(crossprod(T21[, a]))

        #Deflation
        locX1 <- locX1 - T21[, a] %*% t(P1[, a])
        if(!is.null(X3)){
          locX3 <- locX3 - P3[, a] %*% t(T22[, a])
        }
        locX2 <- locX2 - T21[, a] %*% t(P21[, a])
    }

    #Correlations for correlation plot
        R1  <- cor(X1, T21)
        R21 <- t(cor(T21, X2))
        R22 <- t(cor(T22, t(X2)))
        if(!is.null(X3)){
                R3  <- cor(t(X3), T22)
            } else {
                R3 <- matrix(0, 1, npc)
                locX3 <- NULL
            }
  res <- list(call = match.call(),
              maxcomp = a,
              npc = npc,
              R1 = R1,
              R21 = R21,
              R22 = R22,
              R3 = R3,
              T11 = T11,
              T12 = T12,
              T21 = T21,
              T22 = T22,
              T31 = T31,
              T32 = T32,
              P1 = P1,
              P21 = P21,
              P3 = P3,
              var1 = X1,
              var2 = X2,
              var3 = X3,
              var1res = locX1,
              var2res = locX2,
              var3res = locX3)
    class(res) <- c("lplsReg")
    return(res)
}
