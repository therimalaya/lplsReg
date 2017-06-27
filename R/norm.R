#' Vector normalization
#'
#' Normalizes a vector to have length equal to 1.
#'
#'
#' @param vec A numerical vector
#' @return The normalized vector
#' @author Solve Saebo
#' @keywords lpls normalization
#' @examples
#'
#' data(BCdata)
#' x <- norm(BCdata$X[,1])
#' crossprod(x)
#' @export
norm <- function(vec){
  vec / sqrt(drop(crossprod(vec)))
}

#' BCdata Documentation
#' BCdata is the data used in this package
#' @docType data
#' @usage data(BCdata)
#' @name BCdata
#' @docType BCdata
#' @keywords lpls-data BCdata datasets
#' @examples 
#' data(BCdata)
NULL
