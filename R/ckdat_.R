## internal function
## check and format data, extracts n and p for subsequent inputs
##
.ckdat <- function(data) {
  if ("qcdata" %in% class(data)) {
    return(data)
  }

  if (!{
    is.numeric(data) | is.array(data) | is.data.frame(data)
  }) {
    stop("'data' must be a vector, an array, a matrix or a data.frame")
  }


  if (is.array(data) | is.matrix(data)) {
    if (length(dim(data)) > 2) {
      stop("'data' cannot be an array of dimension > 2")
    }
  } else if (is.data.frame(data)) {
    if (!all(sapply(data, is.numeric))) {
      stop("'data' cannot contain non-numeric features")
    } else {
      data <- as.matrix(data)
    }
  } else if (is.numeric(data)) {
    data <- matrix(data, ncol = 1)
  }


  if (any(is.nan(data) | is.infinite(data))) {
    stop("'data' contains NaN and/or Inf values")
  }

  if (any(is.na(data))) {
    stop("'data' contains NA data points. The user is suggested to process NAs with an appropriate method before calling this function.")
  }

  ans <- list(
    data = data,
    unrolled = c(data, recursive = TRUE, use.names = FALSE),
    n = nrow(data),
    p = ncol(data)
  )
  class(ans) <- "qcdata"

  return(ans)
}
