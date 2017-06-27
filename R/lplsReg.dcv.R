#' Regression/classification using lpls with double cross-validation
#'
#' Prediction or classification using lpls and double-cross-validation with
#' potential jackknife variable selection and re-fitting of the lpls model to a
#' reduced variable set based on the inner CV-loop results before prediction of
#' samples in the outer CV-loop.
#'
#' @importFrom stats model.matrix predict
#' @param X1 A response vector or matrix for regression. For classification
#' this should be either a factor or a dummy coded 0/1 matrix with one column
#' per group.
#' @param X2 Predictor matrix of size (n x p).
#' @param X3 Background information matrix of size (m x p)
#' @param npc.sel A vector of component numbers to be tested in the initial
#' LPLS model based on all variables in the inner CV-loop.  Default is 1:5.
#' @param alpha A vector of alpha-values to be tested in the initial LPLS model
#' based on all variables in the inner CV-loop.  Default is a single value 0.
#' See \code{lplsReg} for details on alpha.
#' @param npc.ref A vector of component numbers to be tested in the re-fitted
#' LPLS model based on selected variables in the inner CV-loop.  Default is
#' 1:5.
#' @param testlevel Testlevel for the jackknife testing of the variables.
#' Deafult is 0.05
#' @param dreduce Logical. Should variable selection on the columns of X3
#' (parallel to X2) also be applied to the rows of X3? This is logical only if
#' X3 is a (p x p) matrix expressing some dependency or simlarity between the
#' variables in X2,  hence,  in cases where both the rows and columns of X3
#' relate to the variables of X2.
#' @param colcent Logical vector of length referring to X2 and X3. Should
#' column centering be performed?
#' @param rowcent Logical vector of length referring to X2 and X3. Should row
#' centering be performed?
#' @param grandcent Logical vector of length referring to X2 and X3. Should
#' overall centering be performed?
#' @param outerfolds A list of length \code{k} defining the sample numbers in
#' each fold of k-fold cross-validation (outer segments).  May use
#' \code{balanced.folds} to make the segment list up front of the analysis.
#' @param innernfolds The number of segments to be used in the inner cv-loop.
#' @param err.type The evaluation criterion for prediction/classification
#' performance. Either "rate" (total error rate),  "rmsep" (root mean square
#' error),  or "rmsep2" a modified rmsep where only predictions between 0 and 1
#' contribute to the error.  Predictions outside this range are considered as
#' perfect predictions.
#' @return
#' \item{call}{The function call}
#' \item{total.error }{The total error according to the
#' chosen evaluation criterion (see \code{err.type})}
#' \item{apost }{For classification: Posterior probabilities of class
#' membership for each sample.}
#' \item{trueclass }{For clasification: The true class of each sample}
#' \item{class }{For clasification: The predicted class of each sample}
#' \item{X1pred }{The numerical predicted value for each sample.}
#' \item{cv.npc.sel }{The best number of components for the the outer-loop
#' model chosen from the inner-loop results. This is used for the outer loop
#' predictions if \code{refit = FALSE}.}
#' \item{cv.alpha}{The best value of alpha
#' for the the outer-loop model chosen from the inner-loop results.}
#' \item{cv.npc.ref }{The best number of components for the the outer-loop
#' model chosen from the results of re-fitted inner-loop models to jackknife
#' selected variables. This is used for the outer loop predictions if
#' \code{refit = TRUE}.}
#' \item{varfreq }{The frequency of variables being
#' selected across outer cross-validation segments.}
#' @author Solve Saebo
#' @keywords lpls double-CV regression classification
#' @examples
#'
#' data(BCdata)
#' segs <- balanced.folds(BCdata$Y,  5)
#' fit.dcv <- lplsReg.dcv(factor(BCdata$Y),  BCdata$X,  BCdata$Z,  
#'                        outerfolds  =  segs,  innernfolds  =  5)
#' @export
lplsReg.dcv <-
#function(X1,  X2,  X3,  npc.sel = 1:5,  alpha = 0,  npc.ref = 1:5,  refit = FALSE,  testlevel = 0.05,  dreduce = F,
  function(X1,  X2,  X3,  npc.sel = 1:5,  alpha = 0,  npc.ref = NULL,  testlevel = 0.05,  dreduce = F,
           colcent = c(T, T), rowcent = c(F, F),  grandcent = c(F, F),  outerfolds,  innernfolds = 10,  err.type = "rate"){
  #Performs double crossvalidation using L-PLS and Jackknifing for variable selection
  #Optimal number components for LPLS-model refitting on jackknife variables is found in inner loop,  and
  #used for prediction of left out-samples from the outer loop.
  #outerX2  =  matrix (n x p) of predictors
  #outerX1  =  classindicator (dummy matrix (0/1) or a factor variable)
  #outerX3  =  background information matrix (l x p)
  #npc.sel  =  the set of LPLS components to evaluate
  #alpha  =  level of background inclusion
  #npc.ref  =  the set of LPSL-components to evaluate for refit on selected variables
  #testlevel  =  testlevel for jackknife-testing of variable significance
  #dreduce  =  should variable selection be performed on both rows,  and columns of X3 (only natural if
  #X3 is square and symmetric).
  #outerfolds  =  Specification of the segments of the outer CV-loop. If "LOO" leave-one-out is used
  #innernfolds  =  Number of segments in the inner CV-loop. If "LOO" leave-one-out is used

    nouter  <- dim(X2)[1]
    pouter  <- dim(X2)[2]
    q1      <- dim(X1)[2]
    dcvX2   <- X2
    dcvX1   <- X1
    dcvX3   <- X3
    varfreq <- rep(0, pouter)

    if (as.character(outerfolds)[1] == "LOO") {outerfolds <- as.list(1:nouter)}
    outernfolds <- length(outerfolds)
    classlabs   <- NULL
    trueclass   <- NULL

    #Checking the response matrix whether prediction or classification is to be done

    if (!is.null(q1) && q1 == 1) {
      q1 <- NULL
      X1 <- X1[, 1, drop = TRUE]
    }
    if (is.null(q1)) {
      if (class(X1) == "factor") {
        dcvX1     <- model.matrix(~X1 - 1)
        classlabs <- levels(X1) #colnames(dcvX1)
        trueclass <- apply(dcvX1, 1, function(x){which(x == 1)})
        q1        <- dim(dcvX1)[2]
        mod       <- "Classification"
      } else {
        dcvX1 <- as.matrix(X1)
        q1    <- dim(dcvX1)[2]
        mod   <- "Regression"
        if (err.type == "rate") {
          err.type <- "rmsep"
          warning("Argument err.eval.type has been changed to 'rmsep'")
        }
      }
    } else if (q1 >1) {
      rsum <- apply(X1, 1, sum)
      if (all(rsum == 1) & all(X1 %in% c(0, 1))) {
        mod <- "Classification"
        classlabs <- colnames(X1)
        if (is.null(classlabs)) {
          classlabs <- as.character(1:q1)
        }
        trueclass <- apply(X1,  1,  function(x){which(x == 1)})
      } else {
        mod <- "Regression"
        if (err.type == "rate") {
          err.type <- "rmsep"
          warning("Argument err.eval.type has been changed to 'rmsep'")
        }
      }
      dcvX1 <- X1
    }

    X1pred      <- matrix(0, nouter, q1)
    bestalpha   <- rep(0, outernfolds)
    bestnpc.sel <- rep(0, outernfolds)
    bestnpc.ref <- rep(0, outernfolds)

    for(jj in 1:outernfolds){
      innX1 <- dcvX1[-outerfolds[[jj]], , drop = F]
      innX2 <- dcvX2[-outerfolds[[jj]], , drop = F]
      outX2 <- as.matrix(dcvX2[outerfolds[[jj]], , drop = F])
      if (innernfolds == "LOO") {
        innfolds <- as.list(1:dim(innX1)[1])
      }
      else{
        innfolds <- balanced.folds(innX1[, 1, drop = F], innernfolds)
      }

      innercv <- lplsReg.cv(
        innX1, innX2, dcvX3, npc.sel = npc.sel,
        alphavek = alpha, npc.ref = npc.ref,
        colcent = colcent,  rowcent = rowcent,
        grandcent = grandcent, dreduce = dreduce,
        folds = innfolds,  testlevel = testlevel,
        err.eval.type = err.type
      )

      #Optimal model identification:
      bestind   <- which(innercv$err == min(innercv$err), arr.ind = TRUE)
      nocomp    <- which(bestind[, 3] == 1)
      refitting <- ifelse(length(nocomp) == 0, TRUE, FALSE)
      if (!refitting) {
        bestind <- bestind[nocomp, , drop = F]
        bestind <- bestind[bestind[, 1] == min(bestind[, 1]), , drop = F]
        bestind <- bestind[bestind[, 2] == max(bestind[, 2]), , drop = F]
        bestind <- bestind[1, , drop = F]
      }
      else{
        bestind <- bestind[bestind[, 3] == min(bestind[, 3]), , drop = F]
        bestind <- bestind[bestind[, 2] == max(bestind[, 2]), , drop = F]
        bestind <- bestind[1, , drop = F]
      }

      bestnpc.sel[jj]   <- npc.sel[bestind[1]]
      bestalpha[jj]     <- alpha[bestind[2]]
      bestnpc.ref[jj]   <- ifelse(refitting, npc.ref[(bestind[3]-1)], 0)
      usecomp           <- ifelse(refitting, bestnpc.ref[jj], bestnpc.sel[jj])

      varindex          <- innercv$sigvars[, , bestind[1], bestind[2]]
      varfreq[varindex] <-varfreq[varindex] + 1/outernfolds
      if (dreduce) {
        innX3 <- dcvX3[varindex, varindex, drop = F]
      } else {
        innX3 <- dcvX3[, varindex, drop = F]
      }
      fit <- lplsReg(
        X1 = innX1, X2 = innX2[, varindex], X3 = innX3,
        npc = usecomp,  alpha = 0.5,
        rowcentering = rowcent,  colcentering = colcent,
        grandcentering = grandcent
      )
      predfit <- predict(fit, usecomp, X2new = outX2[, varindex, drop = F])
      X1pred[outerfolds[[jj]], ] <- predfit$pred
      cat("Outer segment", jj, "completed \n")
    }
    classpred <- NULL
    apost     <- NULL
    if (err.type == "rate") {
      apost     <- t(apply(X1pred, 1, function(x){exp(x) / sum(exp(x))}))
      classpred <- apply(apost, 1, function(x){which(x == max(x))})
      err       <- confusion(trueclass, classpred, classlabs)
    }
    else if (err.type == "rmsep") {
      err <- sqrt(mean((X1pred - X1) ^ 2))
    }

    res <- list(
      call        = match.call(),
      total.error = err,
      apost       = apost,
      trueclass   = trueclass,
      class       = classpred,
      X1pred      = X1pred,
      cv.npc.sel  = bestnpc.sel,
      cv.alpha    = bestalpha,
      cv.npc.ref  = bestnpc.ref,
      varfreq     = varfreq
    )
    class(res) <- "lplsReg.dcv"
    return(res)
  }
