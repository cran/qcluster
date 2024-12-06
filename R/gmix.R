#
#  gmix
#
#  Wrapper for ECM in C
#
#  Arguments:
#  - (data.frame/list) data            : data points (N, P)
#  - (integer)           K             : the number of clusters
#  - (char/matrix/fn)    init          : {'kmed', 'kmeans', 'pam'}, matrix, vector or callable (see notes)
#  - (numeric)           erc           : eigen-ratio constraint for gmix
#  - (ingteger)          iter.max      : maximum number of iterations for gmix
#  - (numeric)           tol           : threshold to determine convergence of gmix
#  - (integer)           init.nstart   : number of times kmedian is run to determine the best solution
#  - (numeric)           init.tol      : threshold to determine convergence of kmedian
#  - (logic)             save_cluster  : should point-to-cluster assignment be saved in output?
#  - (logic)             save_params   : should estimated gmix triplets be saved in output?
#  - (logic)             save_taus     : should fuzzy point-to-cluster assignment be saved in output?
#
#  Returns:
#  - (list) ans            : a list containing results of gmix algorithm (10 elmts)
#              $info       : information on gmix convergence and errors (see notes)
#              $iter       : number of iterations run by gmix algo
#              $N          : number of data points
#              $P          : data dimension
#              $K          : number of clusters
#              $loglik    : expected loglikelihood from last iteration of gmix
#              $size       : cluster sizes
#              $cluster    : point-to-cluster assignment
#              $taus       : fuzzy point-to-cluster assignment
#              $params     : list of clusters' triplets
#                     $prop: clusters' mixing proportions (pi's)
#                     $mean: clusters' centers
#                     $cov : clusters' covariances
#
#  Notes:
#  - init: (if results in error the function returns NA)
#     'kmed'            : uses kmedians (kmedians.c)
#     'kmeans' or 'pam' : cluster data with kmeans or pam and transforms the clustering vector
#                         into a set of taus weight (0-1)
#     matrix            : must be a matrix of initial weights of dimensions (N, K)
#     vector            : must be a vector labels (N); converted into matrix of weights
#     function          : must be a function that takes in input data and number of clusters (K)
#                         and returns a matrix of taus weight (as matrix above)
#  - The list in output from C is refined in R

#  - The ans$info value is composed of two informations code + flag:
#    (In R it is converted by .gmix_info function)
#
#    flags (embedded in "info")
#        0 = noflag
#        1 = numerically degenerate tau cannot be prevented
#        2 = successfully enforced the ERC at least once
#        3 = 1 + 2
#
#    info: 3 digits: 1st digit (hundreds) is for code
#      -200: used to signal DSYEV fail (returns NA)  (code=-2 in R)
#      -100: used to signal MALLOC fail (returns NA) (code=-1 in R)
#       10X: no better than initial                  (code= 3 in R)
#       20X: iter.max reached                        (code= 2 in R)
#       30X: convergence within iter.max             (code= 1 in R)
#        (X): where X is from flags.
#
#
gmix <- function(
    data, K = NA, erc = 50, iter.max = 1000, tol = 1e-8,
    init = "kmed", init.nstart = 25, init.iter.max = 30,
    init.tol = tol, save_cluster = TRUE, save_params = TRUE, save_taus = FALSE) {
  # Check data in input
  data <- .ckdat(data)

  N <- data$n
  P <- data$p
  unroldat <- data$unrolled
  data <- data$data

  # Check K
  if (!(is.na(K) || (!is.na(K) && length(K) == 1 && K == floor(K) && (K > 0)))) {
    stop("'K' must be a positive integer or NA")
  }

  # Set up initialization
  ecm_init <- NULL
  if (is(init, "character")) {
    if (init == "kmed") {
      init <- TRUE
      ecm_init <- NaN
    } else if (init == "kmeans") {
      init <- FALSE
      try(silent = TRUE, {
        ecm_init <- kmeans(data, centers = K, iter.max = init.iter.max, nstart = init.nstart)$cluster
      })
    } else if (init == "pam") {
      init <- FALSE
      try(silent = TRUE, {
        ecm_init <- pam(data, k = K, cluster.only = TRUE)
      })
    } else {
      # Try calling function (must have signature function(data,K) -> tau matrix)
      init <- FALSE
      try(silent = TRUE, {
        ecm_init <- do.call(get(init), data, K)
      })
    }
    if (!init) {
      if (!is.null(ecm_init)) {
        ecm_init <- as.numeric(sapply(seq(K), function(x) ecm_init == x))
        init <- FALSE
      } else {
        message("Initialization failed. Returning NA")
        return(NA)
      }
    }
  } else if (is(init, "matrix") | is(init, "data.frame")) {
    if (!is.na(K) && K != ncol(init)) {
      stop(sprintf("number of clusters in 'init' (=%d) and 'K' (=%d) ", ncol(init), K))
    } else {
      ecm_init <- as.numeric(init)
      init <- FALSE
    }
  } else if (is(init, "vector")) {
    if (length(init) != N) {
      stop(sprintf("length(init)=%d; nrow(data)=%d; if 'init' is a vector, it must provide a label for each data point", N, length(init)))
    } else {
      if (is.na(K)) {
        K <- length(unique(init))
      } else if (K != length(unique(init))) {
        stop(sprintf("number of clusters in 'init' (=%d) and 'K' (=%d) ", length(unique(init)), K))
      }
      ecm_init <- as.numeric(sapply(seq(K), function(x) init == x))
      init <- FALSE
    }
  } else if (is(init, "function")) {
    # The function must have signature function(data,K) -> tau matrix
    try(silent = FALSE, {
      ecm_init <- do.call(init, list(data, K))
      init <- FALSE
    })
    if (is.null(ecm_init)) {
      message("Initialization with ", quote(init), " failed. Returning NA")
      return(NA)
    }
  } else {
    stop("Invalid 'init': character, matrix, vector or function are accepted, but ", class(init), " was given.")
  }


  if (is.na(K)) {
    stop("'K' must be given if 'init' is not a \"vector\", \"matrix\", or a \"data.frame\"")
  }
  # Call C ecm
  ans <- .Call(
    C_ECM_C,
    unroldat, N, P, K, init,
    ecm_init, erc, iter.max, tol,
    init.nstart, init.iter.max, init.tol,
    save_cluster, save_params, save_taus
  )


  names(ans)[1:5] <- c("info", "iter", "N", "P", "K")
  ans$info <- .gmix_info(ans$info)

  # Prepare ans
  if (ans$info$code %in% c(1, 2)) {
    names(ans)[6:10] <- c("loglik", "size", "cluster", "posterior", "params")
    # Return loglik
    ans[[6]] <- ans[[6]] * ans[["N"]]
    names(ans$params) <- c("proportion", "mean", "cov")

    # Shift cluster to start from 1
    ans[["cluster"]] <- ans[["cluster"]] + 1

    if (save_params) {
      ans$params$cov <- array(ans$params$cov, dim = c(P, P, K))
    }
  }

  class(ans) <- "mbcfit"
  return(ans)
}
