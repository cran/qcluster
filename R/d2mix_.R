.d2mix <- function(x, y, params) {
  X <- cbind(x, y)
  K <- length(params$proportion)
  n <- length(x)
  prob <- matrix(nrow = n, ncol = K)
  for (k in seq(K)) {
    mu <- params$mean[, k]
    inv_sigma <- solve(params$cov[, , k])
    Delta <- sweep(X, 2, mu, "-")
    smd <- .rowSums(Delta %*% inv_sigma * Delta, length(x), 2)
    prob[, k] <- 1 / (2 * pi) * det(inv_sigma)^0.5 * exp(-smd / 2)
  }
  prob <- sweep(prob, 2, params$proportion, "*")
  prob <- .rowSums(prob, n, K)
  return(prob)
}
