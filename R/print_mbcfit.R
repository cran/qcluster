print.mbcfit <- function(x, ...) {
  code_msg <- switch(as.character(x$info$code),
    "-2" = "Lapack DSYEV failed: not able to compute eigen-decomposition of some cluster covariance matrix",
    "-1" = "Memory allocation error: not enough space for array allocation",
    "1" = "",
    "2" = "'gmix' did not converge (iterations == iter.max)",
    "3" = "EM failed; no better than initial"
  )

  #  flag_msg <- switch(as.character(x$info$flag),
  #    "0" = "",
  #    "1" = "Degenerate posterior-weights could not be prevented",
  #    "2" = "Succesfully Enforced ERC at least once",
  #    "3" = "(Degenerate posterior-weights could not be prevented) AND (Succesfully Enforced ERC at least once)"
  #    )


  if (x$info$code %in% c(-2, -1, 3)) {
    sol_msg <- paste("'gmix' could not find a solution based on the current settings.",
      paste("info: code = ", x$info$code, "; flag = ", x$info$flag, " (See help \"Details\" section.)", sep = ""),
      code_msg,
      sep = "\n"
    )
  } else {
    sol_msg <- paste(paste0("'gmix' discovered K=", x$K, " clusters of size:"), paste(x$size, collapse = " "), sep = "\n")
    sol_msg <- ifelse(code_msg == "", sol_msg, paste(code_msg, sol_msg, sep = "\n"))
  }
  message(sol_msg, "\n")
  message("Available components of the output list:")
  cat(paste(as.character(names(x)), collapse = ", "), end = "\n")
}
