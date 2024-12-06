# For doc see bqs.R
.bqs <- function(data, methodset, B, type, ncores, oob, saveparams) {
  ## Setting null "foreach" running variables for CRAN's checks
  b <- m <- NULL
  #############################################################

  M <- length(methodset)
  dopar <- ifelse((ncores > 1) & ((B > 1) | (M > 1)), TRUE, FALSE)

  # Register parallel backend if needed
  if (dopar) {
    clst <- makeCluster(ncores, type = "PSOCK")
    registerDoParallel(clst)
  } else {
    registerDoSEQ()
  }

  # Export in parallel loops
  parexport <- c()
  parpackgs <- c("qcluster")

  # Add packages require by user_functions
  addpkg <- unique(unlist(lapply(methodset, `[[`, ".packages")))
  parpackgs <- c(parpackgs, addpkg)

  # Add exports required by user_functions
  addexp <- unique(unlist(lapply(methodset, `[[`, ".export")))
  parexport <- c(parexport, addexp)


  # Prepare bootstrap indexes
  if (B == 0) {
    bootidx <- matrix(seq(data$n), ncol = 1)
    B <- 1
  } else {
    bootidx <- replicate(B, sample(seq(data$n), size = data$n, replace = TRUE))
  }

  unrolled_bootidx <- matrix(seq(0, data$p - 1) * data$n, ncol = data$p, nrow = data$n, byrow = TRUE)
  unrolled_bootidx <- apply(bootidx, 2, function(x) c(sweep(unrolled_bootidx, 1, x, "+")))

  # Prepare type switch
  type <- switch(type,
    both = 2L,
    smooth = 1L,
    hard = 0L
  )

  ret <- list(boot_id = bootidx)


  ############ Scoring Utility ##################################
  .nocheck_qscore <- function(unroldat, N, P, params, type) {
    K <- length(params$proportion)
    ans <- .Call(C_SCORE_C, type, unroldat, N, P, K, params)
    names(ans) <- c("hard", "smooth")
    return(ans)
  }
  ###############################################################


  # Dedicated runner for oob=0 (no out-of-bag estimates)
  if (oob == 0) {
    container <- foreach(m = 1:M, .combine = c) %:%
      foreach(b = 1:B, .combine = c, .packages = parpackgs, .export = parexport) %dopar% {
        code <- 1L
        if (saveparams) params <- list()
        try(silent = TRUE, expr = {
          params <- methodset[[m]]$fn(data$data[bootidx[, b], ], only_params = TRUE)
          code <- 0L
        })

        score <- c(NA, NA)
        if (code == 0) {
          code <- 2L
          try(silent = TRUE, expr = {
            score <- .nocheck_qscore(data$unrolled, data$n, data$p, params, type)
            code <- 0L
          })
        }
        if (saveparams) {
          list(c(code, score, use.names = FALSE), params)
        } else {
          c(code, score, use.names = FALSE)
        }
      }
    if (saveparams) {
      paramslist <- container[seq(2, length(container), by = 2)]
      container <- c(container[seq(1, length(container), by = 2)], recursive = TRUE, use.names = FALSE)
      paramslist <- split(paramslist, sort(rep(seq(M), B)))
      names(paramslist) <- names(methodset)
      ret[["params"]] <- paramslist
    }
    container <-
      array(container,
        dim = c(3, B, M),
        dimnames = list(type = c("code", "hard", "smooth"), boot = NULL, method = names(methodset))
      )
    container <- aperm(container, c(3, 1, 2))
    ret[["scores"]] <- container
  }

  # Dedicated runner for oob=1 (add out-of-bag estimates)
  else if (oob == 1) {
    container <- foreach(m = 1:M, .combine = c) %:%
      foreach(b = 1:B, .combine = c, .packages = parpackgs, .export = parexport) %dopar% {
        score <- c(NA, NA)
        oobscore <- c(NA, NA)

        code <- c(1L, 1L)
        if (saveparams) params <- list()
        try(silent = TRUE, expr = {
          params <- methodset[[m]]$fn(data$data[bootidx[, b], ], only_params = TRUE)
          code <- c(0L, 0L)
        })

        score <- c(NA, NA)
        oobscore <- c(NA, NA)
        if (all(code == 0)) {
          code <- c(2L, 2L)

          try(silent = TRUE, expr = {
            score <- .nocheck_qscore(data$unrolled, data$n, data$p, params, type)
            code[1] <- 0L
          })
          try(silent = TRUE, expr = {
            dt <- data$unrolled[-unrolled_bootidx[, b]]
            oob_N <- as.integer(length(dt) / data$p)
            oobscore <- .nocheck_qscore(dt, oob_N, data$p, params, type)
            code[2] <- 0L
          })
        }
        if (saveparams) {
          list(c(code[1], score, code[2], oobscore, use.names = FALSE), params)
        } else {
          c(code[1], score, code[2], oobscore, use.names = FALSE)
        }
      }
    if (saveparams) {
      paramslist <- container[seq(2, length(container), by = 2)]
      container <- c(container[seq(1, length(container), by = 2)], recursive = TRUE, use.names = FALSE)
      paramslist <- split(paramslist, sort(rep(seq(M), B)))
      names(paramslist) <- names(methodset)
      ret[["params"]] <- paramslist
    }
    # Organize and reshape output
    idx <- matrix(seq(B * M * 6), nrow = 3)
    oob_container <- c(container[c(idx[, seq(2, dim(idx)[2], 2)])])
    container <- c(container[c(idx[, seq(1, dim(idx)[2], 2)])])
    oob_container <-
      array(oob_container,
        dim = c(3, B, M),
        dimnames = list(type = c("code", "hard", "smooth"), boot = NULL, method = names(methodset))
      )
    container <- array(container,
      dim = c(3, B, M),
      dimnames = list(type = c("code", "hard", "smooth"), boot = NULL, method = names(methodset))
    )
    oob_container <- aperm(oob_container, c(3, 1, 2))
    container <- aperm(container, c(3, 1, 2))
    ret[["scores"]] <- container
    ret[["oob_scores"]] <- oob_container
  }


  # Dedicated runner for oob=2 (include ONLY out-of-bag estimates)
  else if (oob == 2) {
    container <- foreach(m = 1:M, .combine = c) %:%
      foreach(b = 1:B, .combine = c, .packages = parpackgs, .export = parexport) %dopar% {
        code <- 1L
        if (saveparams) params <- list()
        try(silent = TRUE, expr = {
          params <- methodset[[m]]$fn(data$data[bootidx[, b], ], only_params = TRUE)
          code <- 0L
        })

        score <- c(NA, NA)
        if (code == 0) {
          code <- 2L
          try(silent = TRUE, expr = {
            dt <- data$unrolled[-unrolled_bootidx[, b]]
            oob_N <- as.integer(length(dt) / data$p)
            score <- .nocheck_qscore(dt, oob_N, data$p, params, type)
            code <- 0L
          })
        }
        if (saveparams) {
          list(c(code, score, use.names = FALSE), params)
        } else {
          c(code, score, use.names = FALSE)
        }
      }
    if (saveparams) {
      paramslist <- container[seq(2, length(container), by = 2)]
      container <- c(container[seq(1, length(container), by = 2)], recursive = TRUE, use.names = FALSE)
      paramslist <- split(paramslist, sort(rep(seq(M), B)))
      names(paramslist) <- names(methodset)
      ret[["params"]] <- paramslist
    }
    container <-
      array(container,
        dim = c(3, B, M),
        dimnames = list(type = c("code", "hard", "smooth"), boot = NULL, method = names(methodset))
      )
    container <- aperm(container, c(3, 1, 2))
    ret[["oob_scores"]] <- container
  }

  if (dopar) {
    stopCluster(clst)
  } else {
    stopImplicitCluster()
  }

  return(ret)
}
