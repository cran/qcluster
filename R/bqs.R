#
#  bqs
#
#  Wrapper for .bqs function caller
#
#  Arguments:
#  - (data.frame) data: The dataset
#  - (list/function) methodset: either a single function or a list of clustering functions
#  - (int/array of int) B: the number of bootsrap replicates (see notes)
#  - (character) type: the type of Quadratic Score in output ("hard", *"smooth"*, "both")
#  - (int) oob: specify whether to save out-of-bag bootstrap results, in {0,1,2} (see notes)
#  - (int) ncores: the number of corse used in foreach (default detectCores - 2); -1 set max core
#  - (float) alpha: confidence-level for empirical bootstrap quantiles (both-tails)
#  - (character) rankby:
# 		- lower-quantile: rank methods by descending lower-bootstrap-quantile (see alpha)
# 		- mean: rank methods by descending mean
# 		- NA: do not rank methods
#  - (float) boot_na_share: methods resulting in more than 'B*boot_na_share' NA estimates are excluded from comparison
#  - (bool) savescores: save boostrap estimated scores in output
#  - (bool) saveparams: save boostrap estimated cluster parameters in output
#
#  Returns:
#  - bqs out
#
#  Notes:
#  - B:
#      0  : fit and score on full data (no bootstrap, no parallel)
#      >0 : fit on bootstrap score on full data (uses parallel)
#
#   - oob:
#          FALSE : ignore out-of-bag scores
#           TRUE : additionally compute and return out-of-bag scores
#         "only" : compute and return ONLY out-of-bag scores (not scoring on full data)
#
bqs <- function(data, methodset, B = 10, type = "smooth",
                oob = FALSE, ncores = detectCores() - 2,
                alpha = 0.05, rankby = ifelse(B==0, "mean", "lq"),
                boot_na_share = 0.25, savescores = FALSE, saveparams = FALSE) {
  # Check data
  data <- .ckdat(data)

  # Check methodset (if a function, wrap it in a list to comply to later code)
  methodset <- mbind(methodset)

  # Check B input
  if (!(length(B) == 1)) {
    stop("bqs: 'B' must be an integer >= 0")
  } else {
    if (!(as.integer(B) == B)) {
      stop("bqs: 'B' must be an integer >= 0")
    }
  }

  # Check type input
  if (!(is.character(type) & length(type) == 1) | !(type[1] %in% c("hard", "smooth", "both"))) {
    stop("bqs: 'type' must be in {'hard','smooth', 'both'}")
  }

  # Check oob input
  if (identical(oob, FALSE)) {
    oob <- 0
  } else if (identical(oob, TRUE)) {
    oob <- 1
  } else if (identical(oob, "only")) {
    oob <- 2
  } else {
    stop("bqs: 'oob' must be one of {TRUE, FALSE, 'only'}")
  }
  if ((oob > 0) & B == 0) {
    stop("bqs: out-of-bag scores can't be computed when fitting clustering on full data; try setting 'oob=FALSE' or 'B>0'")
  }

  # Check alpha
  if (length(alpha) != 1) {
    stop("bqs: 'alpha' must be a number in [0,1]")
  } else if (!is.numeric(alpha)) {
    stop("bqs: 'alpha' must be a number in [0,1]")
  } else if ((alpha > 1) | (alpha < 0)) {
    stop("bqs: 'alpha' must be a number in [0,1]")
  }

  # Check rankby
  criteria <- c("lq", "mean", "1se")
  chkrank <- sapply(criteria, identical, rankby)
  if (!is.na(rankby) && sum(chkrank) != 1) {
    stop("'rankby' must be 'NA' or one of {", paste(paste0("'", criteria, "'"), collapse = ", "), "}")
  }

  # Check save all
  if (!identical(savescores, TRUE) & !identical(savescores, FALSE)) {
    stop("bqs: 'savescores' must be either TRUE of FALSE")
  }

  # Check boot_na_share
  if (length(boot_na_share) != 1) {
    stop("bqs: 'boot_na_share' must be a number in [0,1]")
  } else if (!is.numeric(boot_na_share)) {
    stop("bqs: 'boot_na_share' must be a number in [0,1]")
  } else if ((boot_na_share > 1) | (boot_na_share < 0)) {
    stop("bqs: 'boot_na_share' must be a number in [0,1]")
  }


  # Handle number of cores
  maxCores <- detectCores()
  if (ncores == -1) {
    ncores <- maxCores
  } else if (ncores > maxCores) {
    message("'ncores' is higher than available cores (", maxCores, "); set to maximum available")
    ncores <- maxCores
  } else if (ncores < -1) {
    message("'ncores' must be greater than '-1'; set to 1")
    ncores <- 1
  }

  res <- NULL
  try({
    res <- .bqs(data, methodset, B, type, ncores, oob, saveparams)
  })

  if (!is.null(res)) {
    M <- length(methodset)

    process_res <- function(code_score_array) {
      Bmod <- ifelse(B==0, 1, B)
      scores <- code_score_array[, 2, , drop = FALSE]
      errors <- .rowSums(code_score_array[, "code", , drop = FALSE] != 0, M, Bmod)
      missing <- .rowSums(is.na(scores), M, Bmod)
      nobs <- Bmod - missing
      mn <- .rowMeans(scores, M, Bmod, na.rm = TRUE)
      if (B > 0) {
        sterr <- apply(scores, 1, sd, na.rm = TRUE)

        H <- sweep(scores, 1, mn, "-")
        H <- sweep(H, 1, sqrt(nobs), "*")
        lo <- mn - apply(H, 1, quantile, probs = 1 - alpha / 2, names = FALSE, na.rm = TRUE) / sqrt(nobs)
        up <- mn - apply(H, 1, quantile, probs = alpha / 2, names = FALSE, na.rm = TRUE) / sqrt(nobs)
      } else {
        sterr <- NA
        lo <- NA
        up <- NA
      }

      # List of model in output
      out <- data.frame(
        id = seq(M),
        rank = NA,
        mean = mn,
        sterr = sterr,
        lower_qnt = lo,
        upper_qnt = up,
        n_obs = nobs,
        n_missing = missing,
        row.names = names(methodset)
      )

      return(out)
    }

    out <- list()
    if (oob != 2) {
      if (type %in% c("hard", "both")) {
        out$hard <- process_res(res$scores[, c("code", "hard"), , drop = FALSE])
      }
      if (type %in% c("smooth", "both")) {
        out$smooth <- process_res(res$scores[, c("code", "smooth"), , drop = FALSE])
      }
    }
    if (oob != 0) {
      if (type %in% c("hard", "both")) {
        out$oob_hard <- process_res(res$oob_scores[, c("code", "hard"), , drop = FALSE])
      }
      if (type %in% c("smooth", "both")) {
        out$oob_smooth <- process_res(res$oob_scores[, c("code", "smooth"), , drop = FALSE])
      }
    }
  }

  class(out) <- "bqs"

  out$B <- B
  out$data <- data
  out$methodset <- methodset

  out$rankby <- rankby
  if (!is.na(rankby)) {
    out <- bqs_rank(out, rankby, boot_na_share)
  }

  if (savescores) {
    out$raw <- res[names(res)[names(res) != "params"]]
  }
  if (saveparams) {
    out$raw$params <- res$params
  }

  return(out)
}
