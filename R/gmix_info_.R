#
#  .gmix_info
#
#  convert the C info variable in code + flag for output in R
#
#  Arguments:
#  - int info: info value from C_EMC_C
#
#  Returns:
#  - list info:
#       - info$code: -2 DSYEV Lapack routine failed
#                    -1 dynamic memory allocation fail (malloc|calloc)
#                     1 everything ok
#                     2 not reached convergence (iter==itermax)
#                     3 failure: no better than initial
#
#       - info$flag:  1: numerically degenerate taus
#                     2: successfully enforced ERC at least once
#                     3: 1 + 2
#
#  Notes:
#  - In C, the code part is flipped: 100 = 3 in R, 200 = 2 in R, 300 = 1 in R
#
.gmix_info <- function(info) {
  code <- floor(info / 100)
  code <- ifelse(code < 0, code, 4 - code)
  flag <- info %% 100
  info <- list(code=code, flag=flag)
  return(info)
}

