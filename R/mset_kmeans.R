#
# Return a list of kmeans functions
#
mset_kmeans <- function(K = c(1:10), iter.max = 50, nstart = 30,
                      algorithm = "Hartigan-Wong", trace = FALSE) {
  # Create list of methods
  configs <- expand.grid(
    trace = trace,
    algorithm = algorithm,
    nstart = nstart,
    iter.max = iter.max,
    centers = K,
    stringsAsFactors = FALSE
  )

  args <- rev(colnames(configs))
  args[1] <- "K"
  M <- nrow(configs)
  mlist <- list()

  for (m in seq(nrow(configs))) {
    x <- rev(configs[m, ])

    fulnm <- paste0("kmeans:", paste(args, x, sep = "=", collapse = "|"))


    y <- as.list(x)

    callee <- eval(substitute(
      {
        function(data, only_params = FALSE) {
          res <- do.call(kmeans, c(list(x = data), y))
          res[["params"]] <- clust2params(data, res$cluster)
          if (only_params) {
            return(res$params)
          } else {
            return(res)
          }
        }
      },
      list(y = y)
    ))

    mlist[[m]] <- list(fullname = fulnm, callargs = y, fn = callee)
  }

  maxK <- length(K)
  repK <- M / maxK
  idx <- rep(seq(repK), maxK)
  names(mlist) <- paste0("kmeans_K", configs$centers, "_", idx)

  class(mlist) <- "qcmethod"
  return(mlist)
}
