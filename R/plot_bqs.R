#
#  plot.bqs
#
#  Produce a plot for bqs object ranking results
#
#  Arguments:
#  - x: A bqs object containing ranking results
#  - score: A character vector indicating which scores to plot (default: NULL)
#  - perc_scale: Logical indicating whether to use percentage scale (default: FALSE)
#  - top: Number of top models to highlight (default: NULL)
#  - annotate: Logical indicating whether to annotate the plot (default: NULL)
#  - ...: Additional graphical parameters
#
#  Returns:
#  - Produces a plot and returns NULL invisibly
#
#  Notes:
#  - Valid scores include "hard", "smooth", "oob_hard", "oob_smooth"
#  - The input bqs object must be ranked using bqs_rank function
#  - Additional graphical parameters can be passed through ...
#
plot.bqs <- function(x, score = NULL, perc_scale = FALSE, top = NULL, annotate = NULL, ...) {
  xname <- deparse(substitute(x))
  if (is.na(x$rankby)) {
    stop(xname, " is not ranked; Use bqs_rank to rank")
  }

  valid_score <- c("hard", "smooth", "oob_hard", "oob_smooth")
  if (is.null(score)) {
    score <- valid_score[valid_score %in% names(x)]
  } else {
    if (!(is(score, "character") && all(score %in% valid_score))) {
      stop("Invalid 'score'. Valid scores: ", paste(valid_score, collapse = ", "))
    }
  }

  M <- length(x$methodset)
  if (is.null(top)) {
    origtop <- NULL
  } else if (is(top, "numeric") && (length(top) == 1) && (top <= M)) {
    origtop <- top
  } else {
    stop(sprintf(
      "'top' must be a number lower or equal than 'length(%s$methodset)' (=%d)",
      xname, M
    ))
  }

  if (is.null(annotate)) {
    annotate <- ifelse(M <= 30, TRUE, FALSE)
  } else if (!is(annotate, "logical")) {
    stop("'annotate' must be either TRUE or FALSE")
  }

  # Grep other parameters
  args <- list(...)
  if (is.null(args$col)) {
    col <- "gray40"
  } else {
    col <- args$col
  }

  crit <- switch(x$rankby,
    "lq" = "lower_qnt",
    "1se" = "-1se",
    "mean" = "mean"
  )

  # Find how many graphs are needed
  ng <- 0
  for (sc in score) {
    if (!is.null(x[[sc]]) && (sc %in% score)) {
      ng <- ng + 1
    }
  }
  nr <- ceiling(ng / 2)
  nc <- ceiling(ng / nr)

  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mfrow = c(nr, nc), mar = c(3.1, 4.1, 2.1, 1.1), oma = c(1, 1, 2, 1), mgp = c(1.5, 0.5, 0))
  for (sc in score) {
    res <- x[[sc]]
    res <- res[!is.na(res$rank), ]
    res <- res[order(res$rank, decreasing = FALSE), ]
    M <- nrow(res)
    val <- res[, crit]

    if (is.null(origtop)) {
      top <- sum((val / max(val)) > .95)
    } else {
      top <- origtop
    }

    if (perc_scale) {
      val <- val / max(val) * 100
    } else {
      val <- val - min(val) + 1
    }



    # Use rect to plot bars. Compute vertex.
    stdh <- 1.5
    magh <- 3.5
    offs <- 0.7

    # Compute y-vertex with magnification of top elements
    yb <- c(0)
    yt <- c(stdh)

    if (top > 0) {
      topm <- rev(seq(M))[1:top]
    } else {
      topm <- c()
    }

    for (m in seq(2, M)) {
      # y-vertex with manification of top elements
      if (!(m %in% topm)) {
        yb <- c(yb, yt[m - 1] + offs)
        yt <- c(yt, yb[m] + stdh)
      } else {
        yb <- c(yb, yt[m - 1] + offs)
        yt <- c(yt, yb[m] + magh)
      }
    }

    # Set x-vertex
    xl <- numeric(M)
    xr <- rev(val) # The order of yt, yb is inverted (better is in last positions)

    # Plot
    padded_range <- function(vec) {
      PADDING <- 0.04
      rng <- range(vec)
      return(rng + c(-1, 1) * PADDING * rng)
    }

    if (perc_scale) {
      minx <- min(xr)
      if (annotate) {
        if (minx < 0 && minx > -20) {
          xrange <- c(-20, 115)
        } else if (minx < 0 && minx < -20) {
          xrange <- c(minx, 115)
        } else if (minx > 0 && minx < 20) {
          xrange <- c(-15, 115)
        } else {
          xrange <- c(0, 115)
        }
      } else {
        xrange <- padded_range(xr)
      }
    } else {
      xrange <- padded_range(c(0, xr))
    }
    yrange <- range(c(yt, yb))

    plot(
      x = NA, y = NA, xlim = xrange, ylim = yrange,
      bty = "n", yaxt = "n", xaxt = "n", xlab = "", ylab = ""
    )

    rect(xl, yb, xr, yt, border = "black", col = col)
    grid(nx = NULL, ny = NA)

    # Add title
    title(main = paste0(sc, " score"), cex.main = 1)

    # Readjust tick values on x-axis
    if (perc_scale) {
      axis(1, labels = TRUE)
    } else {
      axis(1, at = axTicks(1), labels = round(axTicks(1) + min(res[[crit]]) - 1, 3), tick = TRUE)
    }

    # Find model id  add on the y-axis for top models
    ytk_lbl <- rev(res$id)
    ytk_at <- ((yt + yb) / 2)
    axis(2,
      at = ytk_at, labels = ytk_lbl, las = 2, cex.axis = 0.8, line = -0.5,
      lwd = 0, lwd.ticks = 0, tck = -0.015, gap.axis = 0.1
    )
    title(ylab = "Methodset ID")

    if (annotate) {
      # Add text with the score (outside bar)
      l_txt <- rev(round(res[, crit], 3))
      x_txt <- xr
      x_txt[x_txt < 0] <- 0
      y_txt <- (yt + yb) / 2
      s_txt <- rep(0.70, M)
      f_txt <- rep(2, M)
      if (top > 0) {
        s_txt[topm] <- rep(0.9, top)
        f_txt[-topm] <- 1
      }
      text(x = x_txt, y = y_txt, label = l_txt, cex = s_txt, pos = 4, offset = 0.2, font = f_txt, col = "black")

      # Find the model name and write into rectangle
      l_txt <- rev(row.names(res))
      x_txt <- xl
      p_txt <- rep(4, M)
      c_txt <- rep("white", M)
      # Modify position and coloring for annotations with negative or short bars
      if (perc_scale) {
        x_txt[(xr > 0) & (xr < 20)] <- 0
        x_txt[(xr < 0) & (xr > -20)] <- xr[(xr < 0) & (xr > -20)]
        p_txt[xr < 20] <- 2
        c_txt[xr < 20 & xr > -20] <- "black"
      }
      text(x = x_txt, y = y_txt, label = l_txt, cex = s_txt, pos = p_txt, offset = 0.2, font = f_txt, col = c_txt)
    }
  }

  mtext(text = sprintf(
    "(%s, ranked by = %s)",
    ifelse(x$B == 0, "no bootstrap", sprintf("B = %s", x$B)),
    crit
  ), side = 3, outer = TRUE, cex = 1.2, font = 2)
}
