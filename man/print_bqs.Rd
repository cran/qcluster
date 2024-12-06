\name{print.bqs}
\alias{print.bqs}

\title{Display Information on Bootstrap Quadratic Scores Objects}

\description{
  This function provides a print method for objects of class \code{bqs}, which are produced by the \code{bqs} function. It prints a summary of the bootstrapped quadratic score results for the clustering solutions considered.
}

\usage{
\method{print}{bqs}(x, ...)
}

\arguments{
  \item{x}{
    An object of class \code{bqs}, usually the output of the \code{bqs} function.
  }
  \item{...}{
    Additional arguments passed to or from other methods.
  }
}
\value{No return value, called for side effects}
\details{
  The \code{print.bqs} function provides a print method for objects of class \code{bqs}.

  If clustering solutions in \code{bqs} are not ranked, the printing method displays a message to the user signalling it. Otherwise, the printing method shows a summary of the top-6 ranked solutions, in decreasing order, for any available scoring method (this is determined by the \code{oob} argument used in input to the \code{bqs} function. See Details in \code{\link{bqs}}).

  The summary tables for ranked methods has \code{row.names} set to the method's codename, and shows the following information along the columns:
    \describe{
        \item{\code{id}}{Method's index in the \code{methodset} list (see Details in \code{\link{bqs}}).}
        \item{\code{rank}}{Method's rank according to ranking criterion.}
        \item{\code{mean}}{Method's mean (bootstrap) quadratic score.}
        \item{\code{sterr}}{Method's standard error for the (bootstrap) quadratic score.}
        \item{\code{lower_qnt}}{(Only shown for "mean" and "lq" ranking) Method's lower \code{alpha/2}-level quantile of the bootstrap distribution of the quadratic score (\code{alpha} is given in input to \code{bqs} function).}
        \item{\code{upper_qnt}}{(Only shown for "mean" and "lq" ranking) Method's upper \code{alpha/2}-level quantile of the bootstrap distribution of the quadratic score (\code{alpha} is given in input to \code{bqs} function).}
        \item{\code{-1se}}{(Only shown for "1se" ranking) Method's mean (bootstrap) quadratic score minus 1 standard error.}
        \item{\code{-1se}}{(Only shown for "1se" ranking) Method's mean (bootstrap) quadratic score plus 1 standard error.}
    }
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
res <- bqs(dat, mlist, B = 100, type = "both", rankby=NA, ncores = 1,
           oob = TRUE, savescores = TRUE, saveparams = FALSE)

# Methods are not ranked; only available components are shown
res

# Rank method and show summaries
ranked_res <- bqs_rank(res, rankby = "lq", boot_na_share = 0.25)

ranked_res
}
}