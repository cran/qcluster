# Compute the multivariate normal densidy according to params
# X is a (n,p) matrix with observations along rows
.dpnorm <- function(X, mu = numeric(ncol(X)), cov = diag(ncol(X)), prop = 1) {
  N <- nrow(X)
  P <- ncol(X)

  inv_sigma <- solve(cov)
  Delta <- sweep(X, 2, mu, "-")
  smd <- .rowSums(Delta %*% inv_sigma * Delta, N, P)
  prob <- (2 * pi)^(-P / 2) * det(inv_sigma)^0.5 * exp(-smd / 2)
  return(prob * prop)
}
