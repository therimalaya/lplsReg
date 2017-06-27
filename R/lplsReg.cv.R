#' Two - class classification using lpls with crossvalidation
#'
#' Regression/classification using lpls and cross - validation with potential
#' jackknife variable selection and optional refitting of model to selected
#' variables.
#'
#' @importFrom stats model.matrix predict pt
#' @param X1 A response vector or matrix for regression. For classification
#' this should be either a factor or a dummy coded 0/1 matrix with one column
#' per group.
#' @param X2 Predictor matrix of size (n x p).
#' @param X3 Background information matrix of size (m x p)
#' @param npc.sel A vector of component numbers to be tested in the initial
#' LPLS model based on all variables in the inner CV - loop.  Default is 1:5.
#' @param alphavek A vector of alpha - values to be tested in the initial LPLS
#' model based on all variables in the inner CV - loop.  Default is a single
#' value 0. See \code{lplsReg} for details on alpha.
#' @param npc.ref A vector of component numbers to be tested in the re - fitted
#' LPLS model based on selected variables in the inner CV - loop.  Default is
#' NULL which gives no refitting.
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
#' @param folds A list of length \code{k} defining the sample numbers in each
#' fold of k - fold cross - validation.  May use \code{balanced.folds} to make the
#' segment list up front of the analysis.
#' @param err.eval.type The evaluation criterion for prediction/classification
#' performance. Either "rate" (total error rate),  "rmsep" (root mean square
#' error),  or "rmsep2" a modified rmsep where only predictions between 0 and 1
#' contribute to the error.  Predictions outside this range are considered as
#' perfect predictions.
#' @param cvreport Logical. Should an iteration report be printed on screen
#' during the computations?
#' @return
#' \item{X1hatmat}{An array holding predicted X1 - values for each number of
#' components (initial model and refitted) and alpha values.}
#' \item{folds }{The CV - segments used.}
#' \item{coefs.all}{An array holding all estimated
#' regression coefficients for all components (initial model) and alphavalues.}
#' \item{sdcoef }{The standard deviations of the regressions coefficients.}
#' \item{trueclass }{For clasification:True class of sample}
#' \item{pval }{The p - values from jackknife testing of each regression coefficient for all
#' levels of components and alpha.}
#' \item{apost }{For clasification:The posterior probability of each sample to
#' belong to each class in case of classification.}
#' \item{class }{For clasification:The predicted class of each
#' sample for all levels of components and alpha.}
#' \item{err }{The total error (as defined by argument \code{err.eval.type}
#' for all level of components and alpha.}
#' \item{sigvars}{An array of logicals defining wether a variable is
#' found to be significant or not. Significance is given for all levels of
#' components and alpha, }
#' @author Solve Saebo
#' @keywords lpls cross - validation regression classification
#' @examples
#'
#' data(BCdata)
#' segs <- balanced.folds(BCdata$Y, 5)
#' fit.cv <- lplsReg.cv(factor(BCdata$Y),  BCdata$X,  BCdata$Z,  folds = segs)
#' @export
lplsReg.cv <-
  function(X1,  X2,  X3,  npc.sel = 1:5,
           alphavek = seq(0, 1, by = 0.2), npc.ref = NULL,  testlevel = 0.05,
           dreduce = F,  colcent = c(T, T),  rowcent = c(F, F),
           grandcent = c(F, F),  folds, err.eval.type = "rate",  cvreport = TRUE){

    n         <- dim(X2)[1]
    p         <- dim(X2)[2]
    q1        <- dim(X1)[2]
    cvX2      <- X2
    cvX1      <- X1
    cvX3      <- X3

    nfold     <- length(folds)
    classlabs <- NULL
    trueclass <- NULL

    #Checking the response matrix whether prediction or classification is to be done
    if (!is.null(q1) && q1 == 1){
      q1 <- NULL
      X1 <- X1[, 1, drop = TRUE]
    }
    if (is.null(q1)){
      if (class(X1) == "factor"){
        cvX1      <- model.matrix(~X1 - 1)
        classlabs <- levels(X1)#colnames(cvX1)
        trueclass <- apply(cvX1,  1,  function(x){which(x == 1)})
        q1        <- dim(cvX1)[2]
        mod       <- "Classification"
      } else {
        cvX1 <- as.matrix(X1)
        q1   <- dim(cvX1)[2]
        mod  <- "Regression"
        if (err.eval.type == "rate"){
          err.eval.type <- "rmsep"
          warning("Argument err.eval.type has been changed to 'rmsep'")
        }
      }
    } else if (q1 >1) {
      rsum <- apply(X1, 1, sum)
      if (all(rsum == 1) & all(X1 %in% c(0, 1))){
        mod       <- "Classification"
        classlabs <- colnames(X1)
        trueclass <- apply(X1,  1,  function(x){which(x == 1)})
      } else {
        mod <- "Regression"
        if (err.eval.type == "rate"){
          err.eval.type <- "rmsep"
          warning("Argument err.eval.type has been changed to 'rmsep'")
        }
      }
      cvX1 <- X1
    }

    #    if (!refit){npc.ref <- NULL}
    Npc.sel   <- length(npc.sel)
    Nalpha    <- length(alphavek)
    Npc.ref   <- length(npc.ref)
    X1hatmat  <- array(0, dim = c(n, q1, Npc.sel, Nalpha, (Npc.ref + 1)))

    # classes <- apply(cvX1, 1, function(x){which(x == 1)})
    coefs.all <- array(0, dim = c(p, q1, Npc.sel, Nalpha))
    coefsum   <- array(0, dim = c(p, q1, Npc.sel, Nalpha))
    coefqvsum <- array(0, dim = c(p, q1, Npc.sel, Nalpha))

    #Model - fitting and prediction PPLS
    for(j in 1:nfold){
      calX2  <- cvX2[-folds[[j]], , drop = F]
      calX1  <- cvX1[-folds[[j]], , drop = F]
      testX2 <- as.matrix(cvX2[folds[[j]], , drop = F])

      for(m1 in 1:Nalpha){
        calib <- lplsReg(calX1, calX2, cvX3, max(npc.sel), alphavek[m1], colcentering = colcent,
                         rowcentering = rowcent, grandcentering = grandcent)

        #Fitting model to all samples to collect parameter estimates
        if (j == nfold){
          allcalib <- lplsReg(cvX1, cvX2, cvX3, max(npc.sel), alphavek[m1], colcentering = colcent,
                              rowcentering = rowcent, grandcentering = grandcent)
        }
        for(m2 in 1:Npc.sel){
          fit                               <- predict(calib,  ncomp = npc.sel[m2],  X2new = calX2)
          testfit                           <- predict(calib,  ncomp = npc.sel[m2],  X2new = testX2)
          X1hatmat[folds[[j]], , m2, m1, 1] <- testfit$pred
          coefsum[, , m2, m1]               <- coefsum[, , m2, m1]  +  fit$beta1
          coefqvsum[, , m2, m1]             <- coefqvsum[, , m2, m1]  +  fit$beta1^2
          if (j == nfold){
            allfit <- predict(allcalib,  ncomp = npc.sel[m2],  X2new = cvX2)
            coefs.all[, , m2, m1]<-allfit$beta1
          }
        }
      }
      if (cvreport){cat("Segment", j, "of", nfold, "completed \n")}
    }#j

    sdmat <- array(0, dim = c(p, q1, Npc.sel, Nalpha))
    for(m2 in 1:Npc.sel){
      sdmat[, , m2, ] <- sqrt((nfold - 1)/nfold*(coefqvsum[, , m2, ]-(1/nfold)*(coefsum[, , m2, ]^2)))
    }

    #Finding significant variables using Jackknife
    tobs  <- -abs(coefs.all)/sdmat
    pval  <- 2*pt(tobs, df = (nfold - 1))
    issig <- pval<testlevel

    #Refitting and internal crossvalidation on selected variables
    if (!is.null(npc.ref)){
      cat("Internal crossvalidation using JK - variables only \n")
      for(j in 1:nfold){
        for(m1 in 1:Nalpha){
          for(m2 in 1:Npc.sel){
            varind <- which(apply(issig[, , m2, m1, drop = FALSE], 1, any))
            calX1  <- cvX1[-folds[[j]], , drop = F]
            calX2  <- cvX2[-folds[[j]], varind, drop = F]
            testX2 <- as.matrix(cvX2[folds[[j]], varind, drop = F])
            if (dreduce){calX3 <- cvX3[varind, varind, drop = F]}else(calX3 <- cvX3[, varind, drop = F])
            redcalib <- lplsReg(calX1, calX2, calX3, max(npc.ref), alphavek[m1], colcentering = colcent,
                                rowcentering = rowcent, grandcentering = grandcent)
            for(m3 in 1:Npc.ref){
              redtestfit <- predict(redcalib,  ncomp = npc.ref[m3],  X2new = testX2)
              X1hatmat[folds[[j]], , m2, m1, (m3 + 1)]<-redtestfit$pred
            }#m3
          }#m2
        }#m1
        cat("Segment", j, "of", nfold, "completed \n")
      }#j
    }#end refit
    apost     <- NULL
    classpred <- NULL
    if (mod == "Classification"){
      apost     <- array(0,  dim = dim(X1hatmat))
      classpred <- array(0,  dim = c(n, Npc.sel, Nalpha, (Npc.ref + 1)))
    }
    errmat <- array(0, dim = c(Npc.sel, Nalpha, (Npc.ref + 1)))
    for(m1 in 1:Nalpha){
      for(m2 in 1:Npc.sel){
        for(m3 in 1:(Npc.ref + 1)){
          if (err.eval.type == "rate"){
            apost[, , m2, m1, m3]   <- t(apply(X1hatmat[, , m2, m1, m3, drop = F], 1, function(x){exp(x)/sum(exp(x))}))
            classpred[, m2, m1, m3] <- apply(apost[, , m2, m1, m3, drop = F], 1, function(x){which(x == max(x))})
            conf                    <- confusion(trueclass, classpred[, m2, m1, m3], classlabs)
            errmat[m2, m1, m3]      <- conf$total.err
          }
          else if (err.eval.type == "rmsep"){
            errmat[m2, m1, m3] <- sqrt(mean((X1hatmat[, , m2, m1, m3] - cvX1)^2))
          }
        }
      }
    }
    dimnames(errmat) <- list(paste("Comp.sel", npc.sel), paste("alpha", alphavek), paste("Comp.refit", c(0, npc.ref)))
    dimnames(issig)  <- list(NULL, NULL, paste("Comp.sel", npc.sel), paste("alpha", alphavek))
    sigvars          <- issig[, 1, , , drop = F]
    res              <- list(
      X1hatmat = X1hatmat,
      folds = folds,
      coefs.all = coefs.all,
      sdcoef = sdmat,
      trueclass = trueclass,
      pval = pval,
      apost = apost,
      class = classpred,
      err = errmat,
      sigvars = sigvars
    )
    res
  }
