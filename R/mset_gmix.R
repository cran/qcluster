#
# Return a list of gmix functions
#
mset_gmix <- function(K = seq(10), init = "kmed", erc = c(1, 50, 1000),
                     iter.max = 1000, tol = 1e-8,
                     init.nstart = 25, init.iter.max = 30, init.tol = tol) {

  # Handle multiple init types
  if (class(init)[1] == "character") {
    list_init <- as.list(init)
    names(list_init) <- init
  } else if (class(init)[1] %in% c("matrix", "data.frame")) {
    # In case init is a single matrix of weights
    list_init <- list(user_weights = as.matrix(init))
    init <- "user_weights"
  } else if (class(init)[1] == "function") {
    # In case init is a single function
    list_init <- list(user_function = init)
    init <- "user_function"
  } else if (class(init)[1] == "list") {
    # Handle lists of characters, matrix and functions
    nms_init <- character(length(init))
    list_init <- list()
    fn_cnt <- 0
    mt_cnt <- 0
    for (i in seq_along(init)) {
      cls <- class(init[[i]])[1]
      if (cls == "character") {
        nms_init[[i]] <- init[[i]]
        list_init[[i]] <- init[[i]]
      } else if (cls == "function") {
        nms_init[[i]] <- ifelse(fn_cnt > 0, paste0("user_function_", fn_cnt), "user_function")
        list_init[[i]] <- init[[i]]
        fn_cnt <- fn_cnt + 1
      } else if (cls %in% c("matrix", "data.frame")) {
        nms_init[[i]] <- ifelse(mt_cnt > 0, paste0("user_weights_", mt_cnt), "user_weights")
        list_init[[i]] <- init[[i]]
        mt_cnt <- mt_cnt + 1
      }
    }
    names(list_init) <- nms_init
    init <- nms_init
  }


  # Expand arguments into configurations
  configs <- expand.grid(
    init.tol = init.tol,
    init.iter.max = init.iter.max,
    init.nstart = init.nstart,
    tol = tol,
    iter.max = iter.max,
    erc = erc,
    init = init,
    K = K,
    stringsAsFactors = FALSE
  )


  args <- rev(colnames(configs))
  M <- nrow(configs)
  mlist <- list()

  for (m in seq(M)) {
    x <- rev(configs[m, ])

    fulnm <- paste0("gmix:", paste(args, x, sep = "=", collapse = "|"))

    y <- as.list(x)
    y$init <- list_init[[y$init]]

    callee <- eval(substitute(
      {
        function(data, save_cluster = TRUE, save_params = TRUE, save_taus = TRUE, only_params = FALSE) {
          if (only_params) {
            save_cluster <- FALSE
            save_taus <- FALSE
          }
          res <- do.call(gmix, c(
            list(
              data = data, save_cluster = save_cluster,
              save_params = save_params, save_taus = save_taus
            ),
            y
          ))
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
  names(mlist) <- paste0("gmix_K", configs$K, "_", idx)

  class(mlist) <- "qcmethod"
  return(mlist)
}
