#
#  bqs_select
#
#  Function to select solutions from a bqs object based on ranking and type.
#
#  Arguments:
#  - (bqs) bqs_sol: An object of class 'bqs' containing the clustering solution.
#  - (integer) rank: The rank of the solution to select (default 1).
#  - (character) type: The type of Quadratic Score to use; one of {"hard", "smooth", "oob_hard", "oob_smooth"} (default "smooth").
#  - (character) rankby: Criteria to rank solutions by; one of {"lq", "mean", "1se"} or NA to use existing ranking.
#  - (float) boot_na_share: Threshold for excluding methods based on the share of NA estimates in bootstrapping (default 0.25).
#
#  Returns:
#  - (list) A list of solutions that match the given rank and type criteria.
#
#  Notes:
#  - The bqs object must be pre-ranked using the `bqs_rank` function if rankby criteria are provided.
#  - Function halts with an error if bqs_sol is not of class 'bqs', or if invalid 'rank', 'type', or 'rankby' values are provided.
#
bqs_select <- function(bqs_sol, rank = 1, type = "smooth", rankby = NA, boot_na_share = 0.25) {

  bname <- deparse(substitute(bqs_sol))
  if (!inherits(bqs_sol, "bqs")) {
    stop(bname, ' is not of class "bqs"')
  }

  # Check rank
  M <- length(bqs_sol$methodset)
  if (is.numeric(rank) && length(rank) == 1) {
    if (rank < 0) {
      stop("'rank' must be greater than 1")
    }
  } else {
    stop("'rank' must be an integer greater than 1")
  }

  # Check type
  types <- c("hard", "smooth", "oob_hard", "oob_smooth")
  chktype <- sapply(types, identical, type)
  if (sum(chktype) != 1) {
    stop("'type' must be one of {", paste(paste0("'", types, "'"), collapse = ", "), "}")
  }

  # Check rankby
  if (!is.na(rankby)) {
    criteria <- c("lq", "mean", "1se")
    chkrank <- sapply(criteria, identical, rankby)
    if (sum(chkrank) != 1) {
      stop("'rankby' must be one of {", paste(paste0("'", criteria, "'"), collapse = ", "), "}")
    }
    bqs_sol <- bqs_rank(bqs_sol, rankby, boot_na_share)
  } else if (is.na(bqs_sol$rankby)) {
    stop("Solutions not ranked. Use `bqs_rank` to rank.")
  }

  # Check solution is available
  if (!type %in% names(bqs_sol)) {
    stop(type, ' score is not available in ', bname)
  }

  # Determine available ranking to extract (mr)
  sdat <- bqs_sol[[type]]
  maxrank <- min(M, max(sdat$rank, na.rm = TRUE))
  if (maxrank < rank) {
    message(sprintf("Worst-ranked solution for %s scoring is rank %d; returning this as 'rank'=%d is not available",
                    type, maxrank, rank))
    mr <- maxrank
  } else {
    mr <- rank
  }

  # Find all solutions with the given rank
  sol_id <- na.omit(sdat[sdat$rank == mr, 'id'])
  if (length(sol_id) > 0) {
    sol_list <- list()
    message(sprintf("%s score... found %d rank-%d solution(s).", type, length(sol_id), mr), end = "\n")
    # Optionally keep track of rank and rankby method
    # sol_list[["rank"]] <- mr
    # sol_list[["rankby"]] <- bqs_sol$rankby
    for (i in sol_id) {
      nm <- names(bqs_sol$methodset)[i]
      sol_list[[nm]] <- "estimation failed"
      try({
        message(sprintf("\t...estimating solution: %s (methodsed id = %d)", nm, i), end = "\n")
        sol_list[[nm]] <- bqs_sol$methodset[[i]]$fn(bqs_sol$data$data)
      }, silent = TRUE)
    }
    return(sol_list)
  } else {
    message(sprintf("%s score... ! No rank-%d solution found", type, mr))
    return(NULL)
  }
}