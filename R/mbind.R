#
# Check and prepare a list of clustering methods
#
mbind <- function(...) {
  elements <- list(...)
  mlist <- list()
  counter <- 0
  for (el in elements) {
    if (is.function(el)) {
      counter <- counter + 1
      nm <- paste0("user_function_", counter)
      callee <- eval(substitute(function(data, only_params = FALSE, ...) {
        res <- f(data)
        if (only_params) {
          return(res$params)
        } else {
          return(res)
        }
      }), list(f = el))
      mlist[[nm]] <- list(fullname = NULL, callargs = list(f = substitute(el)), fn = callee)
    } else if (identical("qcmethod", class(el))) {
      mlist <- c(mlist, el)
    } else if (identical("list", class(el))) {
      if (!all(sapply(el, is.function))) {
        stop("List input argument must be a list of functions")
      } else {
        for (f in el) {
          counter <- counter + 1
          nm <- paste0("user_function_", counter)
          callee <- eval(substitute(function(data, only_params = FALSE, ...) {
            res <- fn(data)
            if (only_params) {
              return(res$params)
            } else {
              return(res)
            }
          }), list(fn = f))
          mlist[[nm]] <- list(fullname = NULL, callargs = list(f = substitute(f)), fn = callee)
        }
      }
    } else {
      stop("Input argument must be either functions, list of functions, or output from mk_* functions")
    }
  }

  class(mlist) <- "qcmethod"
  return(mlist)
}
