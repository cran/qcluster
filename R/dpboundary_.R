# Compute the classification of points in R^p according to params
# X is a (n,p) matrix with observations along rows
.dpclassify <- function(X, params) {
  prop <- params$prop
  mean <- params$mean
  cov <- params$cov

  K <- length(prop)
  Z <- matrix(nrow = nrow(X), ncol = K)
  for (k in seq(K)) {
    Z[, k] <- .dpnorm(X, mean[, k], cov[, , k], prop[k])
  }
  Z <- apply(Z, 1, which.max)
  return(Z)
}
