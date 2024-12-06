#
# Return a list of pam functions
# NOTE: user can input only parameters that actually modify the algorithm;
#       flags for controlling output are passed to each produced function in mlist
#
mset_pam <- function(K = seq(10),
                     metric = "euclidean",
                     medoids = if (is.numeric(nstart)) "random",
                     nstart = if (variant == "faster") 1 else NA,
                     stand = FALSE,
                     do.swap = TRUE,
                     variant = "original",
                     pamonce = FALSE) {
  # Produce argument for expand.grid in 'expandlist'
  expandlist <- list()
  # variant and pamonce can't be specified together
  # handle user input as in PAM
  if (!missing(variant) & !missing(pamonce)) {
    stop("Set either 'variant' or 'pamonce', but not both")
  } else if (!missing(pamonce)) {
    expandlist[["pamonce"]] <- pamonce
  } else {
    expandlist[["variant"]] <- variant
  }

  expandlist[["do.swap"]] <- do.swap
  expandlist[["stand"]] <- stand

  if (!missing(nstart)) {
    expandlist[["nstart"]] <- nstart
  }

  list_medoids <- NULL
  if (!missing(medoids)) {
    if (is(medoids, "character")) {
      expandlist[["medoids"]] <- medoids
      list_medoids <- list(medoids)
      names(list_medoids) <- names(medoids)
    } else if (is(medoids, "integer") || is(medoids, "numeric")) {
      message("Setting 'K' argument on the base of 'medoids'")
      K <- NULL
      expandlist[["medoids"]] <- "user_labels"
      list_medoids <- list(user_labels = as.integer(medoids))
    } else if (is(medoids, "list")) {
      message("Setting 'K' argument on the base of 'medoids'")
      K <- NULL
      list_medoids <- list()
      nms_medoids <- c()
      lb_cnt <- 0
      for (i in seq_along(medoids)) {
        if (is(medoids[[i]], "character")) {
          nms_medoids <- c(nms_medoids, medoids[[i]])
          list_medoids[[i]] <- medoids[[i]]
        } else if (is(medoids[[i]], "integer") || is(medoids[[i]], "numeric")) {
          nm <- ifelse(lb_cnt > 0, paste0("user_labels_", lb_cnt), "user_labels")
          nms_medoids <- c(nms_medoids, nm)
          list_medoids[[i]] <- as.integer(medoids[[i]])
          lb_cnt <- lb_cnt + 1
        }
      }
      names(list_medoids) <- nms_medoids
      expandlist[["medoids"]] <- nms_medoids
    }
  }

  expandlist[["metric"]] <- metric
  expandlist[["k"]] <- K

  # Make configuration matrix
  configs <- expand.grid(expandlist, stringsAsFactors = FALSE)

  if (!is.null(list_medoids)) {
    configs[["k"]] <- as.numeric(sapply(configs$medoids, function(x) {
      length(list_medoids[[x]])
    }))
  }

  # Make list of models
  args <- rev(colnames(configs))
  args[1] <- "K"
  M <- nrow(configs)
  mlist <- list()

  for (m in seq(nrow(configs))) {
    x <- rev(configs[m, ])

    fulnm <- paste0("pam:", paste(args, x, sep = "=", collapse = "|"))

    y <- as.list(x)
    if (!missing(medoids)) {
      y$medoids <- list_medoids[[y$medoids]]
    }

    callee <- eval(substitute(
      {
        function(data, diss = inherits(data, "dist"), cluster.only = FALSE,
                 keep.data = FALSE, keep.diss = FALSE, only_params = FALSE) {
          if (only_params | cluster.only) {
            cluster.only <- TRUE
          }
          res <- do.call(pam, c(
            list(
              x = data, diss = diss, cluster.only = cluster.only,
              keep.data = keep.data, keep.diss = keep.diss
            ),
            y
          ))
          if (cluster.only) {
            res <- clust2params(data, res)
          } else {
            res[["params"]] <- clust2params(data, res$clustering)
          }
          return(res)
        }
      },
      list(y = y)
    ))

    mlist[[m]] <- list(fullname = fulnm, callargs = y, fn = callee)
  }

  maxK <- length(unique(configs[["k"]]))
  repK <- M / maxK
  idx <- rep(seq(repK), maxK)
  names(mlist) <- paste0("pam_K", configs$k, "_", idx)

  class(mlist) <- "qcmethod"
  return(mlist)
}
