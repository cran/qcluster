# Plot data, optionally adding clustering information
plot_clustering <- function(data,  subset = NULL,
                            cluster = NULL, params = NULL,
                            what = c("clustering", "contour", "boundary"),
                            col_cl = NULL, pch_cl = NULL) {

  # Set main features
  dt <- .ckdat(data)
  data <- dt$data
  P <- dt$p

  # Check inputs
  if (is.null(cluster) && is.null(params)) {
    stop("At least one of 'cluster' or 'params' must be given")
  }

  K <- NA
  if (!is.null(cluster)) {
    if (nrow(data) != length(cluster)) {
      stop(sprintf("'cluster' (obs %d) must be of the same length as data (%d)", length(cluster), nrow(data)))
    }
    K <- length(unique(cluster))
  }

  if (!is.null(params)) {
    .ckclp(params)
    if (is.na(K)) {
      K <- length(params$proportion)
    } else {
      if (K != length(params$proportion)) {
        stop(sprintf("The number of groups from 'cluster' argument (%d) does not match that of 'params' (%d)",
                     K, length(params$proportion)))
      }
    }
    if (P != nrow(params$mean)){
      stop(sprintf("Data dimension (%d) does not match that implied by 'params' (%d)",
                   P, nrow(params$mean)))
    }
  }

  # set pch by cluster
  if (is.null(pch_cl)) {
    pch_cl <- c(1:9, letters, LETTERS)[seq(K)]
  }

  # set colors by cluster
  if (is.null(col_cl)) {
    col_cl <- seq(K)
  }

  if (("clustering" %in% what) && is.null(cluster)) {
    message("'clustering' won't be plotted as 'cluster' was not given")
  }

  if ((("contour" %in% what) || ("boundary" %in% what))) {
    if (is.null(params)) {
      message("'contour' and/or 'boundary' won't be plotted as 'params' was not given")
    }
  }

  # Check subset parameter
  # this is used to focus on two data feature only
  if (!is.null(subset)) {
    if (is(subset, "numeric")) {
      if (length(subset) == P) {
        message(sprintf("'subset' ignored as data is %d-dimensional", P))
      } else {
        data <- data[, subset, drop = FALSE]
        P <- length(subset)
        if (!is.null(params)) {
          params$proportion <- params$proportion
          params$mean <- params$mean[subset, , drop = FALSE]
          params$cov <- params$cov[subset, subset, , drop = FALSE]
        }
      }
    } else {
      stop("'subset' must be a vector indexing columns of 'data'")
    }
  }


  padded_range <- function(vec) {
    PADDING <- 0
    rng <- range(vec)
    rng <- rng + c(-1, 1) * PADDING * rng
    return(rng)
  }

  base_plot <- function(xdata, ydata, xlim = NULL, ylim = NULL) {
    if (is.null(xlim)) {
      xrange <- padded_range(xdata)
    } else {
      xrange <- xlim
    }
    if (is.null(ylim)) {
      yrange <- padded_range(ydata)
    } else {
      yrange <- ylim
    }
    plot(
      x = NA, y = NA, xlim = xrange, ylim = yrange, xlab = "", ylab = "",
      axes = TRUE, frame.plot = TRUE, xaxt = "n", yaxt = "n"
    )
  }

  plot_boundary <- function(params) {
    axlims <- par("usr")
    x <- seq(axlims[1], axlims[2], length.out = 300)
    y <- seq(axlims[3], axlims[4], length.out = 300)
    X <- as.matrix(expand.grid(x, y))
    Z <- matrix(.dpclassify(X, params), nrow = length(x), ncol = length(y))
    image(x, y, Z, col = adjustcolor(col_cl, alpha.f = 0.1), useRaster = TRUE, add = TRUE)
  }

  # PLOTTING
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  if (P == 1) {
    yfake <- -0.02

    # GREP PROPER YLIM if contour is required
    if (("contour" %in% what) && !is.null(params)) {
      x <- range(data)
      x <- seq(x[1], x[2], length.out = 1000)
      Z <- matrix(nrow = length(x), ncol = K)
      for (k in seq(K)) {
        Z[, k] <- dnorm(x, mean = params$mean[, k], sd = sqrt(params$cov[, , k])) * params$proportion[k]
      }
      yl <- range(c(yfake, Z))
    } else {
      yl <- NULL
    }

    # BASE EMPTY PLOT
    data <- cbind(data, yfake)
    par(mar = c(5.1, 2.1, 4.1, 2.1))
    base_plot(data[, 1], data[, 2], ylim = yl)
    title(xlab = colnames(data)[1])
    axis(1, labels = TRUE)
    axlims <- par("usr")

    # PLOT CONTOUR
    if (("contour" %in% what) && !is.null(params)) {
      axis(2, labels = TRUE)
      for (k in seq(K)) {
        lines(x, Z[, k], col = col_cl[k])
      }
    }

    # PLOT BOUNDARY
    if (("boundary" %in% what) && !is.null(params)) {
      # Find points separating regions
      x <- seq(axlims[1], axlims[2], length.out = 1000)
      Z <- matrix(nrow = length(x), ncol = K)
      for (k in seq(K)) {
        Z[, k] <- dnorm(x, mean = params$mean[, k], sd = sqrt(params$cov[, , k])) * params$proportion[k]
      }
      Z <- apply(Z, 1, which.max)
      chgpt <- Z[-1] - Z[-length(Z)]
      chgpt <- x[chgpt != 0]
      xl <- c(axlims[1], chgpt)
      xr <- c(chgpt, axlims[2])
      rect(
        xleft = xl, ybottom = rep(axlims[3], length(col_cl)),
        xright = xr, ytop = rep(axlims[4], length(col_cl)),
        border = NA, col = adjustcolor(col_cl, alpha.f = 0.1)
      )
    }

    # PLOT DATA
    if (("clustering" %in% what) && !is.null(cluster)) {
      for (k in seq(K)) {
        points(data[cluster == k, 1], data[cluster == k, 2], col = col_cl[k], pch = pch_cl[k])
      }
    } else {
      points(data[, 1], data[, 2], pch = 19)
    }
  } else if (P == 2) {
    # BASE EMPTY PLOT
    base_plot(data[, 1], data[, 2])
    title(xlab = colnames(data)[1], ylab = colnames(data)[2])
    axis(1, labels = TRUE)
    axis(2, labels = TRUE)

    # PLOT CONTOUR
    if (("contour" %in% what) && !is.null(params)) {
      for (k in seq(K)) {
        .d2ellipse(params$proportion[k], params$mean[, k], params$cov[, , k], ell_col = col_cl[k])
      }
    }

    # PLOT BOUNDARY
    if (("boundary" %in% what) && !is.null(params)) {
      plot_boundary(params)
    }


    # PLOT DATA
    if (("clustering" %in% what) && !is.null(cluster)) {
      for (k in seq(K)) {
        points(data[cluster == k, 1], data[cluster == k, 2], col = col_cl[k], pch = pch_cl[k])
      }
    } else {
      points(data[, 1], data[, 2], pch = 19)
    }
  } else if (P > 2) {
    par(mfcol = c(P, P), mar = c(.5, .5, .5, .5), oma = c(3, 3, 3, 3))
    for (i in seq(P)) {
      for (j in seq(P)) {
        base_plot(data[, i], data[, j])
        # ADD AXES TICKS ALL AROUND
        if (i == 1) {
          axis(2, labels = TRUE)
        } else if (i == P) {
          axis(4, labels = TRUE)
        }
        if (j == 1) {
          axis(3, labels = TRUE)
        } else if (j == P) {
          axis(1, labels = TRUE)
        }

        # PLOT i,j graph
        if (i == j) {
          axlims <- par("usr")
          text(x = mean(axlims[1:2]), y = mean(axlims[3:4]), colnames(data)[i], cex = 1.5)
        } else {
          subparams <- list()
          subparams$proportion <- params$proportion
          subparams$mean <- params$mean[c(i, j),]
          subparams$cov <- params$cov[c(i, j), c(i, j),]
          # PLOT CONTOUR
          if (("contour" %in% what) && !is.null(params)) {
            for (k in seq(K)) {
              .d2ellipse(subparams$proportion[k], subparams$mean[, k], subparams$cov[, , k], ell_col = col_cl[k])
            }
          }

          # PLOT BOUNDARY
          if (("boundary" %in% what) && !is.null(params)) {
            plot_boundary(subparams)
          }


          # PLOT DATA
          if (("clustering" %in% what) && !is.null(cluster)) {
            for (k in seq(K)) {
              points(data[cluster == k, i], data[cluster == k, j], col = col_cl[k], pch = pch_cl[k])
            }
          } else {
            points(data[, i], data[, j], pch = 19)
          }
        }
      }
    }
  }
}
