plot.mbcfit <- function(x, data = NULL, subset = NULL,
                        what = c("clustering", "contour"),
                        col_cl = NULL, pch_cl = NULL, ...) {

  xname <- deparse(substitute(x))
  if (!inherits(x, "mbcfit")) {
    stop(xname, ' is not of class "mbcfit"')
  }

  if (!x$info$code %in% c(1, 2)) {
    print(x)
    stop(xname, " does not contain enough information to produce a new plot")
  } else {
    if (is.null(data)) {
      data <- as.data.frame(t(cbind(
        x$params$mean + 6 * apply(x$params$cov, 3, diag),
        x$params$mean - 6 * apply(x$params$cov, 3, diag)
      )))
      colnames(data) <- paste0("Var. ", seq(ncol(data)))
      pch_cl <- rep(NA, x$K)
      cluster <- c(seq(x$K), seq(x$K))
      what <- c("clustering", what)
    } else {
      cluster <- x$cluster
    }
    plot_clustering(data = data, cluster = cluster, params = x$params, pch_cl = pch_cl, col_cl = col_cl, subset = subset, what = what)
  }
}
