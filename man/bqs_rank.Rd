\name{bqs_rank}
\alias{bqs_rank}

\title{Ranking Clusters Quadratic Scores Estimated Via Boostrap}

\description{Ranks the scores of clusters methods estimated via boostrap}


\usage{
   bqs_rank(bqsol, rankby = "lq", boot_na_share = 0.25)
   }

\arguments{
  \item{bqsol}{
    an object of class \code{bqs} obtained from \code{\link{bqs}}.
    }
  \item{rankby}{
    character string specifying how the scored solution are ranked.
    Possible values are \code{{"lq", "mean", "1se"}}.
    With \code{="lq"} (default), the solutions are ranked by maximizing the estimated
    lower limit of the of the \code{1-alpha} bootstrap confidence
    intervalfor  the expected score.
    With \code{="mean"}, the solutions are ranked by maximizing the estimated
    expected score. 
    With \code{="1se"}, the solutions are ranked by maximizing the estimated
    lower limit of the confidence interval for the expected score
    whose semi-length is equal to a \emph{standard error}.
    The expected score's \emph{standard error} is approximated using the
    boostrap distribution.
  }
  \item{boot_na_share}{
    a numeric value in \code{(0,1)}.
    During the boostrapping a method's score is set to \code{NA} if the
    underlying comptutation runs into errors. 
    Methods resulting in more than \code{B * boot_na_share} errors are
    excluded from the comparison
  }
}


\value{
  An S3 object of class \code{bqs}. Output components are those of
  \code{\link{bqs}}. See \emph{Value} in \code{\link{bqs}}.
}

   

\references{
  Coraggio, Luca and Pietro Coretto (2023). Selecting the number of
  clusters, clustering models, and algorithms. A unifying approach based
  on the  quadratic discriminant score.
  \emph{Journal of Multivariate Analysis}, Vol. 196(105181), 1-20.
  doi: \doi{10.1016/j.jmva.2023.105181}
}



\seealso{
  \code{\link{bqs}}
}


\examples{
# load data
data("banknote")
dat <- banknote[-1]

## set up methods
## see also help('mset_user')  and related functions 
KM   <- mset_kmeans(K = 3)
GMIX <- mset_gmix(K=3, erc=c(1,100)) 

# combine tuned methods
mlist <- mbind(KM, GMIX)

# perform bootstrap
# change B and ncores to a much larger value in real problems 
res <- bqs(dat, mlist, B = 3, rankby="lq", ncores=1)
res

   
# now change ranking criterion
res2 <- bqs_rank(res, rankby="mean")
res2
}


