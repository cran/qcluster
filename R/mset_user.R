#
# Return a list of functions
# - fname is the function name as character
# - .packages is the name of the packages the function depends on (see ?foreach)
# - .export is the name of additional exports for foreach (see ?foreach)
# - The function name must be defined and callable from the main env
# - The function MUST return a list containing triplets parameters called "params"
# - MUST take as first argument the data in input
#
mset_user <- function(fname, .packages = NULL, .export = NULL, ...) {
  # Check fname
  if (!is(fname, "character") || length(fname) > 1) {
    stop("mset_user: can only pass a single functions via its name (as character)")
  } else {
    if (!is(get(fname), "function")) {
      stop("input 'fname' must be the name of a callable function")
    }
  }
  # Check packages
  if (!is.null(.packages) && !is(.packages, "character")) {
    stop("mset_user: '.packages' must be a (vector) of character(s)")
  }
  # Check export
  if (!is.null(.export) && !is(.export, "character")) {
    stop("mset_user: '.export' must be a (vector) of character(s)")
  }

  # Create list of methods
  if (length(list(...)) > 0) {
    configs <- expand.grid(..., stringsAsFactors = FALSE)
    args <- rev(colnames(configs))
    M <- nrow(configs)
    mlist <- list()
    for (m in seq(M)) {
      x <- rev(configs[m, ])

      fulnm <- paste0(fname, ":", paste(args, x, sep = "=", collapse = "|"))


      y <- as.list(x)
      callee <- eval(substitute(
        {
          function(data, only_params = FALSE) {
            res <- do.call(func, c(list(data), y))
            if (only_params) {
              return(res$params)
            } else {
              return(res)
            }
          }
        },
        list(func = get(fname), y = y)
      ))

      mlist[[m]] <- list(fullname = fulnm, callargs = y, fn = callee)
      if (length(.packages) > 0) {
        mlist[[m]][[".packages"]] <- .packages
      }
      if (length(.export) > 0) {
        mlist[[m]][[".export"]] <- .export
      }
    }

    names(mlist) <- paste0(fname, "_", seq(M))
  } else {
    mlist <- list()
    fulnm <- fname
    callee <- eval(substitute(
      {
        function(data, only_params = FALSE) {
          res <- do.call(func, list(data))
          if (only_params) {
            return(res$params)
          } else {
            return(res)
          }
        }
      },
      list(func = get(fname))
    ))

    mlist[[1]] <- list(fullname = fulnm, callargs = NULL, fn = callee)
    if (length(.packages) > 0) {
      mlist[[1]][[".packages"]] <- .packages
    }
    if (length(.export) > 0) {
      mlist[[1]][[".export"]] <- .export
    }

    names(mlist) <- fname
  }

  class(mlist) <- "qcmethod"
  return(mlist)
}
