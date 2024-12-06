#
#  bqs_rank
#
#  (re)rank a bqs solution
#
#  Arguments:
#  - bqs bqsol: a solution obtained from bqs function
#  - character rankby: a character indicating the criterion to maximize
#                      "lq" - lower quantile
#                      "mean" - mean
#                      "1se" - mean - 1standarerror
#
#  Returns:
#  - bqsol: ranked solution
#
bqs_rank <- function(bqsol, rankby = "lq", boot_na_share = 0.25) {
  # Check bqsol
  if (!is(bqsol, "bqs")) {
    stop("'bqsol' must be an object obtained from 'bqs' function")
  }

  # Check rankby
  criteria <- c("lq", "mean", "1se")
  chkrank <- sapply(criteria, identical, rankby)
  if (sum(chkrank) != 1) {
    stop("'rankby' must be one of {", paste(paste0("'", criteria, "'"), collapse = ", "), "}")
  }
  if (bqsol$B == 0 && rankby != "mean") {
    stop("only 'rankby=mean' is available when bootstrap was not used (B=0)")
  }

  rank_var <- switch(rankby,
    "lq" = "lower_qnt",
    "mean" = "mean",
    "1se" = "sterr",
    NA
  )

  # Check boot_na_share
  if (length(boot_na_share) != 1 || !is.numeric(boot_na_share) ||
    boot_na_share > 1 || boot_na_share < 0) {
    stop("'boot_na_share' must be a number in [0,1]")
  }

  methodset <- bqsol$methodset
  data <- bqsol$data

  # Rank methods using quadratic scores
  for (type in c("hard", "smooth", "oob_hard", "oob_smooth")) {
    if (!is.null(bqsol[[type]])) {
      dt <- bqsol[[type]]

      criterion <- switch(rank_var,
        "sterr" = dt[, "mean"] - dt[, "sterr"],
        dt[, rank_var]
      )

      # Remove frequently-failing from comparison
      idx <- which(dt$n_missing / (dt$n_missing + dt$n_obs) > boot_na_share)
      criterion[idx] <- NA

      dt$rank <- rank(-criterion, na.last = TRUE, ties.method = "min")
      dt$rank[idx] <- NA

      # sort
      dt <- dt[order(dt$rank),]

      bqsol[[type]] <- dt

      if (!is.na(rankby)) {
        try(
          {
            best <- row.names(dt)[dt$rank==1][1]
            best_mtd <- methodset[[best]]$fn(data$data)
            bqsol[[paste0("best_", type)]] <- list(method_name = best, solution = best_mtd)
          },
          silent = TRUE
        )
      }
    }
  }
  bqsol$rankby <- rankby
  return(bqsol)
}
