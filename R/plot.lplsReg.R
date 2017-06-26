#' Correlation loadings plot for lpls regression model
#'
#' This plot function produces a so-called correlation loadings plot. The
#' correlation loadings are scaled versions of the regular loadings (and scores
#' for X2) which make them informative in an overlayed plot. The plot enables
#' better interpretation of the common covariance patterns in the three data
#' matrices.
#'
#' @importFrom car ellipse
#' @importFrom grDevices dev.new
#' @importFrom graphics abline arrows identify matplot plot points text
#' @param x A model x as returned from \code{lplsReg}
#' @param comps a vector of length 2 indicating which components to be plotted,
#' Default is the two first components.
#' @param doplot A logical vector of length 3 indicating wether the correlation
#' loadings for matrices X1, X2 and X3 should be plotted, respectively.
#' Default is TRUe for all.
#' @param level A numerical vector of length 3 where each element take the
#' value 1 (plot symbol only) or 2 (plot variable labels) for matrix X1, X2 and
#' X3, respectively.
#' @param arrow A numerical vector of length 3 where each element take the
#' value 0 (no arrow) or 1 (arrow from origin).
#' @param xlim Limits for the x-axis
#' @param ylim Limits for the y-axis
#' @param samplecol The color used for the symbols for the "correlation scores"
#' of X2.
#' @param pathcol The color used for the correlation loading arrows (if
#' arrow = 1) for X3.
#' @author Solve S<c3><a6>b<c3><b8>
#' @keywords correlation-loadings biplot
#' @examples
#'
#' data(BCdata)
#' fit <- lplsReg(BCdata$Y, BCdata$X, BCdata$Z, npc = 10)
#' plot(fit)
#' @export
plot.lplsReg <- function(x, comps = c(1, 2), doplot = c(TRUE, TRUE, TRUE),
                         level = c(2, 2, 2), arrow = c(1, 0, 1), xlim = c(-1, 1),
                         ylim = c(-1, 1), samplecol = 4, pathcol = 2) {
  plottype <- c("p","n")
  plot(
    x    = xlim,
    y    = ylim,
    type = "n",
    xlab = paste("Comp", comps[1]),
    ylab = paste("Comp", comps[2]),
    main = paste("Correlation loading plot")
  )
  ellipse(c(0, 0), matrix(c(1, 0, 0, 1), 2, 2), radius = 1, lty = 1, lwd = 1, col = 1, center.pch = F)
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)

  if (doplot[2]) {
    points(x$R21[, comps[1]], x$R21[, comps[2]], type = plottype[level[2]], pch = 20, col = "grey70", cex = 2)
    points(x$R22[, comps[1]], x$R22[, comps[2]], type = plottype[level[2]], pch = 20, col = 4, cex = 2)
    if (level[2] == 2) {
      text(x$R21[, comps[1]], x$R21[, comps[2]], labels = colnames(x$var2), cex = 0.7, col = "grey70")
      text(x$R22[, comps[1]], x$R22[, comps[2]], labels = rownames(x$var2), cex = 0.7, col = samplecol)
    }
    if (arrow[2] == 1) {
      arrows(0, 0, x$R21[, comps[1]], x$R21[, comps[2]], lwd = 2, col = "grey70", length  =  0.1)
      arrows(0, 0, x$R22[, comps[1]], x$R22[, comps[2]], lwd = 2, col = 4, length  =  0.1)
    }
  }

  if (doplot[1]) {
    points(x$R1[, comps[1]], x$R1[, comps[2]], type = plottype[level[1]], pch = 20, col = 3, cex = 2)
    if (level[1] == 2) {
      text(x$R1[, comps[1]], x$R1[, comps[2]], labels = colnames(x$var1), cex = 0.7, col = 3)
    }
    if (arrow[1] == 1) {arrows(0, 0, x$R1[, comps[1]], x$R1[, comps[2]], col = 3, length  =  0.1)}
  }

  if (doplot[3]) {
    points(x$R3[, comps[1]], x$R3[, comps[2]], type = plottype[level[3]], pch = 20, col = 2, cex = 2)
    if (level[3] == 2) {
      text(x$R3[, comps[1]], x$R3[, comps[2]], labels = rownames(x$var3), cex = 0.7, col = pathcol)
    }
    if (arrow[3] == 1) {arrows(0, 0, x$R3[, comps[1]], x$R3[, comps[2]], col = pathcol, length  =  0.1)}
  }
}
