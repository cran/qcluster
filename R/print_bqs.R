print.bqs <- function(x, ...) {
  if (is.na(x$rankby)) {
    cat("Solutions not ranked. Use `bqs_rank` to rank.", end = "\n")
  } else {
    for (type in c("hard", "smooth", "oob_hard", "oob_smooth")) {
      if (!is.null(x[[type]])) {
        # Prepare prompt string
        bootstr <- ifelse(x$B == 0, "no bootstrap", sprintf("B=%d", x$B))
        if (x$B > 0) {
          rankbystr <- switch(x$rankby,
            "lq" = ", rankby='lq'",
            "1se" = ", rankby='1se'",
            "mean" = ", rankby='mean'"
          )
        } else {
          rankbystr <- ""
        }
        line1 <- sprintf("%s score (%s%s):", type, bootstr, rankbystr)
        line2 <- paste(rep("-", nchar(line1)), collapse = "")
        cat(line1, line2, sep = "\n")

        # Select at most 6 entries and hide/rename
        # columns according to B and rankby
        dt <- x[[type]]
        dt <- dt[!is.na(dt$rank), ]
        dt <- dt[order(dt$rank), ]
        showup <- min(6, nrow(dt))

        if (x$B == 0) {
          dt <- dt[seq(showup), c("id", "rank", "mean")]
          colnames(dt) <- c("id", "rank", "score")
        } else {
          if (x$rankby == "1se") {
            dt <- dt[seq(showup), c("id", "rank", "mean", "sterr")]
            dt[["-1se"]] <- dt$mean - dt$sterr
            dt[["+1se"]] <- dt$mean + dt$sterr
            dt <- dt[, c("id", "rank", "mean", "sterr", "-1se", "+1se")]
          } else {
            dt <- dt[seq(showup), c(
              "id", "rank", "mean", "sterr", "lower_qnt",
              "upper_qnt"
            )]
          }
        }
        print(format(dt, digits = 4))
        cat("\n\n")
      }
    }
  }
  cat("\nAvailable components:\n", names(x), "\n\n")
}
