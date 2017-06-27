#' Matrix centering
#'
#' A function which performs row-,  column-,  double- or overall-centering of a
#' matrix
#'
#' If both row- and colcenter are \code{TRUE} the matrix is double centered.
#' This is done by first subtracting the pre-computed row-means,  then the
#' pre-computed column means and finally,  the grand mean is added back in.
#'
#' @param Mat A numerical matrix.
#' @param rowcenter Logical. If \code{TRUE} the matrix is row-centered
#' @param colcenter Logical. If \code{TRUE} the matrix is columns centered
#' @param grandcenter Logical. If \code{TRUE} the matrix is only centered by a
#' subtraction by the overall (grand) mean of all matrix elements. Only applies
#' if both rowcenter and colcenter are \code{FALSE}.
#' @return A centered matrix
#' @author Solve Saebo
#' @keywords lpls centering
#' @examples
#'
#' data(BCdata)
#' Xc <- centering(BCdata$X,  FALSE,  TRUE,  FALSE)
#' @export
centering <- function(Mat, rowcenter = FALSE,  colcenter = FALSE,  grandcenter = FALSE){
  #Function for centering of matrices

  dims <- dim(Mat)

  if (rowcenter & colcenter) {
    rowm   <- apply(Mat, 1, mean)
    colm   <- apply(Mat, 2, mean)
    grandm <- mean(Mat)
    Mat    <- Mat - t(matrix(rep(1, dims[2]), ncol = 1) %*% rowm) -
      matrix(rep(1, dims[1]), ncol = 1) %*% colm  +
      matrix(grandm, nrow = dims[1], ncol = dims[2])
  }
  else if (rowcenter) {
    Mat    <- t(scale(t(Mat), scale = F))
    rowm   <- attr(Mat, "scaled:center")
    colm   <- rep(0, dims[2])
    grandm <- 0
  }
  else if (colcenter) {
    Mat    <- scale(Mat, scale = F)
    colm   <- attr(Mat, "scaled:center")
    rowm   <- rep(0, dims[1])
    grandm <- 0
  }
  else if (grandcenter) {
    rowm   <- rep(0, dims[1])
    colm   <- rep(0, dims[2])
    grandm <- mean(Mat)
    Mat    <- Mat - matrix(grandm, nrow = dims[1], ncol = dims[2])
  }
  else{
    rowm   <- rep(0, dims[1])
    colm   <- rep(0, dims[2])
    grandm <- 0
  }
  attr(Mat, "rowm")   <- rowm
  attr(Mat, "colm")   <- colm
  attr(Mat, "grandm") <- grandm
  Mat
}
