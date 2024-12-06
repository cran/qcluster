# Draw 50% CI and 95% CI normal ellipses of 'cov', centered on mean

.d2ellipse <- function(prop, mean, cov, add = TRUE, ell_col = "black") {
  INTERVAL50 <- 1.17741 # sqrt of quantile of chisquare df=2 at alpha level
  INTERVAL95 <- 2.447747 # or equiv: sqrt(-2*log(1-alpha))
  # Spectral decomposition of covariance matrix
  eig <- eigen(cov)
  eigvec <- eig$vectors
  eigval <- eig$values
  # Unit circle
  t <- seq(0, 2 * pi, length.out = 100)
  X <- cbind(cos(t), sin(t))
  # First: scale by sqrt of eigenvalues
  # (sqrt, because we square to compute the covariance of scaled circle)
  X <- sweep(X, 2, sqrt(eigval), "*")
  # Second: rotate by eigenvectors
  T <- X %*% t(eigvec)
  # Third: scale to stay within 50 and 95 CI (i.e. multiplying by
  # the relative quantile)
  T50 <- T * INTERVAL50
  T95 <- T * INTERVAL95
  T50 <- sweep(T50, 2, mean, "+")
  T95 <- sweep(T95, 2, mean, "+")

  if (!add) {
    plot(mean[1], mean[2],
      xlim = range(T95[, 1], finite = TRUE),
      ylim = range(T95[, 2], finite = TRUE),
      type = "n", asp = 1,
      xlab = "", ylab = ""
    )
  }

  # Plot contours
  lines(T50[, 1], T50[, 2], col = ell_col, lty = "dashed", lwd=1.5)
  lines(T95[, 1], T95[, 2], col = ell_col, lty = "dotted", lwd=1.5)
  # Add text. Location is obtained finding a point on the unit circle
  # and then using the (covariance) affine transformation and scaling by INTERVAL* to
  # obtain the quantile
  t50at <- (c(cos(pi / 2), sin(pi / 2)) * sqrt(eigval)) %*% t(eigvec) * INTERVAL50 + mean
  t95at <- (c(cos(pi / 2), sin(pi / 2)) * sqrt(eigval)) %*% t(eigvec) * INTERVAL95 + mean
  # Grep current background color
  bg <- par("bg")
  bg <- ifelse(bg == "transparent", "white", bg)
  # Add text (grep the probability at the quantile level, scaled by proportion)
  lab50 <- round(.dpnorm(cbind(t50at[1], t50at[2]), mean, cov, prop), 3)
  lab95 <- round(.dpnorm(cbind(t95at[1], t95at[2]), mean, cov, prop), 3)
  # Grep text width and height (once plot is called)
  tw <- max(strwidth(lab50, cex = 0.75), strwidth(lab95, cex = 0.75))
  th <- max(strheight(lab50, cex = 0.75), strheight(lab95, cex = 0.75))
  # Draw rect to draw text background
  rect(
    xleft = t50at[1] - tw / 2, ybottom = t50at[2] - th / 2,
    xright = t50at[1] + tw / 2, ytop = t50at[2] + th / 2,
    col = adjustcolor(bg, alpha.f = 1), border = NA
  )
  rect(
    xleft = t95at[1] - tw / 2, ybottom = t95at[2] - th / 2,
    xright = t95at[1] + tw / 2, ytop = t95at[2] + th / 2,
    col = adjustcolor(bg, alpha.f = 1), border = NA
  )
  # Plot label text
  text(x = t50at[1], y = t50at[2], cex = 0.75, labels = lab50, col = adjustcolor(ell_col, alpha.f = 0.75))
  text(x = t95at[1], y = t95at[2], cex = 0.75, labels = lab95, col = adjustcolor(ell_col, alpha.f = 0.75))
}
