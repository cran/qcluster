\name{bqs_select}
\alias{bqs_select}

\title{Select Ranked Cluster Solutions by Quadratic Score}

\description{
  Select solutions from a \code{bqs} object based on specified rank and type of score.
}

\usage{
   bqs_select(bqs_sol, rank = 1, type = "smooth", rankby = NA, boot_na_share = 0.25)
}

\arguments{
  \item{bqs_sol}{
    An object of class \code{bqs} containing the clustering solutions to be selected.
  }
  \item{rank}{
    An integer \code{>0} specifying the rank of the solution to select. Default is \code{1}.
  }
  \item{type}{
    A character string specifying the type of Quadratic Score. Possible values are \code{"hard"}, \code{"smooth"}, \code{"oob_hard"}, \code{"oob_smooth"}. Default is \code{"smooth"}.
  }
  \item{rankby}{
    A character string specifying the criteria used to rank solutions in \code{bqs}. Possible values are \code{"lq"}, \code{"mean"}, \code{"1se"}, or \code{NA} (default). See \emph{Details}.
  }
  \item{boot_na_share}{
    A numeric value between (0, 1). Clustering solutions in \code{bqs_sol} with a share of \code{NA} bootstrap estimates are excluded from ranking. Default is \code{0.25}.
  }
}

\value{
  A named list of all clustering solutions achieving a \code{type} score of rank \code{rank} when ranked  according to \code{rankby} criterion, or \code{NULL} if no such solution is available in the \code{bqs_sol} object. List names correspond to methods' names, and each named entry contains the corresponding clustering method in \code{bqs_sol$methodlist} fit on \code{bqs_sol$data}.
}

\details{
Even if the \code{bqs_sol} object is not pre-ranked, the user may specify a ranking criterion to rank clustering solutions dynamically using the \code{rankby} argument; this does not influence the \code{bqs_sol} ranking. In these instances, the user can also specify \code{boot_na_share} as in \code{\link{bqs_rank}} to exclude solutions based on the proportion of unsuccessful bootstrap estimations. If \code{rankby=NA}, the \code{bqs_sol} must be pre-ranked.
}

\seealso{
  \code{\link{bqs}}, \code{\link{bqs_rank}}
}

\examples{
\donttest{
# Load data and set seet
set.seed(123)
data("banknote")
dat <- banknote[-1]

# set up kmeans, see help('mset_kmeans')
KM    <- mset_kmeans(K = 2:5)

# set up Gaussian model-based clustering via gmix()
GMIX  <- mset_gmix(K=2:5, erc=c(1, 50 , 100))

# combine tuned methods
mlist <- mbind(KM, GMIX)

# perform bootstrap
# se 'ncores' to the number of available physical cores
res <- bqs(dat, mlist, B = 20, type = "both", rankby=NA, ncores = 1,
           oob = TRUE, savescores = TRUE, saveparams = FALSE)

# Methods are not ranked; this will raise an error
try(bqs_select(res, rank = 1))

# Rank method dynamically
ranked_res <- bqs_select(res, rank = 2, rankby = "lq", boot_na_share = 0.25)
names(ranked_res)
}
}

