#
#  score
#
#  Wrapper for score_hard,score_smooth,score_hard_cluster,score_smooth_cluster in C
#
#  Arguments:
#  - (data.frame/matrix) data : data points (N, P)
#  - (list)                cl : a list containing cluster solution triplets
#  - (string)            type : "smooth", "hard" or "both"
#
#  Returns:
#  - (vector) scores          : a vector containing either smooth score, hard score, or both
#

qscore <- function(data, params, type = "both") {
  # Check type input
  type <- switch(type,
    both = 2L,
    smooth = 1L,
    hard = 0L,
    2L
  )

  # Check data in input
  data <- .ckdat(data)

  # Check params in input
  .ckclp(params)

  N <- data$n
  P <- data$p
  unroldat <- data$unrolled

  # Call C score
  K <- length(params$proportion)
  scores <- .Call(C_SCORE_C, type, unroldat, N, P, K, params)
  names(scores) <- c("hard", "smooth")
  return(scores)
}
