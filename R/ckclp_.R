## internal function
## check the consistency of cluster parameters
##
.ckclp <- function(clp) {
  if (!is.list(clp)) {
    stop("Cluster parameters should be provided as a list containing the objects '$proportion', '$mean' and 'cov'")
  }

  if (!{
    "proportion" %in% names(clp) & "mean" %in% names(clp) & "cov" %in% names(clp)
  }) {
    stop("Cluster parameters should be provided as a list containing the objects '$proportion', '$mean' and 'cov'")
  }

  if (!is.numeric(clp$proportion)) {
    stop("The object '$proportion' must be a numeric vector")
  }

  if (!is.array(clp$mean)) {
    stop("The object '$mean' must be a matrix/array where '$mean[,k]' is the mean vector of the k-th cluster")
  }

  if (!is.array(clp$cov)) {
    stop("The object '$cov' must be a 3-dimensional array where '$cov[ , ,k]' is the covariance matrix of the k-th cluster")
  }

  k1 <- length(clp$proportion)
  dm <- dim(clp$mean)
  dc <- dim(clp$cov)

  if (!{
    length(unique(c(k1, dm[2], dc[3]))) == 1
  }) {
    stop("Number of clusters (K) is inconsistent across the objects '$proportion', '$mean' and '$cov'")
  }

  if (!{
    length(unique(c(dm[1], dc[1], dc[2]))) == 1
  }) {
    stop("Data dimension (p) is inconsistent across the objects '$mean' and '$cov'")
  }

  if (abs(sum(clp$proportion) - 1) > sqrt(.Machine$double.eps)) {
    stop("Cluster proportions in '$proportion' should sum up to 1")
  }
}
