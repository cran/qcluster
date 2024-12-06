
\name{mbind}

\alias{mbind}

\title{Combines Methods Settings}

\description{
  The function combines functions containing clustering methods setups 
  built using \code{\link{mset_user}} and related functions.
}

\usage{
  mbind(...)
}

\arguments{
  \item{...}{
    one or more object of class \code{qcmethod} obtained from \code{\link{mset_user}} and related functions
  }
}



\value{
  An S3 object of class \code{'qcmethod'}. Each element of the list
  represents a competing method containing the following objects
  
  \item{fullname}{
    a string identifying the setup.
  }
  \item{callargs}{
    a list with arguments that are passed to the base function.
  }
  \item{fn}{
    the function implementing the specified setting. This \code{fn}
    function can be executed on the data set.
    It has two arguments: \code{data} and \code{only_params}.
    \code{data} is a data matrix or data.frame
    \code{only_params} is logical.
    If \code{only_params==FALSE} (default), \code{fn} will return the
    object returned by the \code{fname}.
    If \code{only_params==TRUE} (default) \code{fn} will return only cluster 
    parameters (proportions, mean, and cov, see \link{clust2params}).
  }
}





\references{
  Coraggio, Luca and Pietro Coretto (2023). Selecting the number of
  clusters, clustering models, and algorithms. A unifying approach based
  on the  quadratic discriminant score.
  \emph{Journal of Multivariate Analysis}, Vol. 196(105181), 1-20.
  doi: \doi{10.1016/j.jmva.2023.105181}
}




\seealso{
  \code{\link{mset_user}}, \code{\link{mset_gmix}},  \code{\link{mset_kmeans}}, \code{\link{mset_pam}}
}







\examples{
# load data
data("banknote")
dat  <- banknote[-1]

# generate  kmeans setups 
A <- mset_kmeans(K=c(2,3))

# generate gmix  setups 
B <- mset_gmix(K=c(2,3))

# combine setups
M <- mbind(A, B)

# get the PAM setting with K=3
m <- M[[4]]
m

# cluster data with M[[3]]
fit <- m$fn(dat)
fit
}

