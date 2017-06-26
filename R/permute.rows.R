#' Internal function (copied from previous version of the package pamr)
#'
#' Function for column permutation within rows of a matrix
#'
#' @importFrom stats runif
#' @param x A numerical matrix
#' @return A matrix with permuted columns within rows.
#' @keywords permutation lpls
#' @examples
#'
#' data(BCdata)
#' permute.rows(BCdata$X[1:10,1:5])
#' @export
permute.rows <- function (x)
{
  dd <- dim(x)
  n  <- dd[1]
  p  <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}
