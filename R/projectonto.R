#' Internal function for projections
#'
#' Projects the (centred) columns of a matrix onto a (centered) vector and
#' returns the slope-coefficient for each column.
#'
#'
#' @param A A centered matrix whose columns will be projected onto
#' @param b this centered vector
#' @return A vector of slope coefficients.
#' @author Solve S<c3><a6>b<c3><b8>
#' @keywords lpls projection
#' @examples
#'
#'
#' data(BCdata)
#' Y <- BCdata$Y
#' X <- BCdata$X
#' #Project the columns of X onto Y.
#' tt <- projectonto(X,Y)
#' @export
projectonto <- function(A,b){
  t(A) %*% b / drop(crossprod(b))
}
