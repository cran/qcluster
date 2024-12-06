#
#  clust2params
#
#  Derive cluster tripltes from data and clusters
#
#  Arguments:
#  - (qcdata) data: data in input
#  - (int array) cluster: integer array specifying point-to-cluster assignment for data
#
#  Returns:
#  - (list) params: a list of cluster's proportion, mean, cov (conformable for scoring)
#
#  Note: For "historical" reasons 'params' are named 'triplets' in C code; (naming as in paper)
#
clust2params <- function(data, cluster){

  # Check data
  data <- .ckdat(data)

  # Check cluster
  if (!is.vector(cluster)){
    stop("'cluster' must be a vector of integers")
  } else if (length(cluster)!=data$n){
    stop("'cluster' length and number of points in 'data' do not match")
  } else {
    cl <- as.integer(factor(cluster)) - 1L
    K <- max(cl) + 1
  }

  params <- .Call(C_TRIPLETS_C, data$unrolled, data$n, data$p, K, cl)
  if (!is.null(params)){
    names(params) <- c("proportion", "mean", "cov")
    params$cov <- array(params$cov, dim = c(data$p, data$p, K))
  }
  return(params)

}
