
\name{mset_gmix}

\alias{mset_gmix}

\title{Generates Methods Settings for Gaussian Mixture Model-Based Clustering}


\description{
  The function generates a software abstraction of a list of clustering
  models implemented through a set of tuned methods and algorithms. 
  In particular, it  generates a list of\code{\link{gmix}} -type functions each
  combining model tuning parameters and other algorithmic settings.
  The generated functions are ready to be called on the data set.
}

\usage{
mset_gmix(
   K = seq(10),
   init = "kmed",
   erc = c(1, 50, 1000),
   iter.max = 1000,
   tol = 1e-8,
   init.nstart = 25, 
   init.iter.max = 30,
   init.tol = tol)
}


\arguments{
  \item{K}{
    a vector/list, specifies the number of clusters.  
  }
  \item{init}{
    a vector, contains the settings of the \code{init} parameter of \code{\link{gmix}}.
  }
  \item{erc}{
    a vector/list, contains the settings of the \code{erc} parameter of \code{\link{gmix}}.
  }
  \item{iter.max}{
    a integer vector, contains the settings of the \code{iter.max} parameter of \code{\link{gmix}}. 
  }
  \item{tol}{
    a vector/list, contains the settings of the \code{tol} parameter of \code{\link{gmix}}. 
  }
  \item{init.nstart}{
    a integer vector, contains the settings of the \code{init.start}
    parameter of \code{\link{gmix}}. 
  }
  \item{init.iter.max}{
    a integer vector, contains the settings of the \code{init.iter.max} parameter of \code{\link{gmix}}. 
  }   
  \item{init.tol}{
    a vector/list, contains the settings of the \code{init.tol} parameter of \code{\link{gmix}}. 
  }   
}


\details{
  The function produces functions implementing competing clustering methods
  based on several Gaussian Mixture models specifications.  
  The function produces functions for fitting competing Gaussian Mixture
  model-based clustering methods settings.
  This is a specialized version of the more general function
  \code{\link{mset_user}}. 
  In particular, it produces a list of \code{\link{gmix}} functions each
  corresponding to a specific setup in terms of both model
  hyper-parameters (\emph{e.g.} the number of clusters, the eigenvalue ratio
  constraint, \emph{etc.}) and algorithm's control parameters
  (\emph{e.g.}  the type of initialization, maximum number of iteration,
  \emph{ etc.}). See \code{\link{gmix}} for a detailed description of
  the   role of each argument and their data types.

}



\value{
  An S3 object of class \code{'qcmethod'}. Each element of the list
  represents a competing method containing the following objects
  
  \item{fullname}{
    a string identifying the setup.
  }
  \item{callargs}{
    a list with  \code{\link{gmix}} function arguments.
  }
  \item{fn}{
    the function implementing the specified setting. This \code{fn}
    function can be executed on the data set.
    It has two arguments: \code{data} and \code{only_params}.
    \code{data} is a data matrix or data.frame
    \code{only_params} is logical.
    If \code{only_params==FALSE} (default), \code{fn} will return the
    object returned by \code{\link{gmix}}.
    If \code{only_params==TRUE} (default) \code{fn} will return only
    cluster parameters (proportions, mean, and cov, see \link{clust2params}.
  }
}







\references{
  Coraggio, Luca, and Pietro Coretto (2023).
  Selecting the Number of Clusters, Clustering Models, and Algorithms.
  A Unifying Approach Based on the Quadratic Discriminant Score.
  \emph{Journal of Multivariate Analysis}, Vol. 196(105181), pp. 1-20,
  \doi{10.1016/j.jmva.2023.105181}
}



\seealso{
  \code{\link{gmix}}, \code{\link{mset_user}}, \code{\link{bqs}}
}







\examples{
# 'gmix' settings combining number of clusters K={3,4} and eigenvalue 
# ratio constraints {1,10} 
A <- mset_gmix(K = c(2,3), erc = c(1,10))
   
# select setup 1: K=2, erc = 1, init =" kmed"
ma1 <- A[[1]]
print(ma1)

# fit M[[1]] on banknote data
data("banknote")
dat  <- banknote[-1]
fit1 <- ma1$fn(dat)   
fit1

# if only cluster parameters are needed
fit1b <- ma1$fn(dat, only_params = TRUE)   
fit1b

   
# include a custom initialization, see also help('gmix')
compute_init <- function(data, K){
  cl  <- kmeans(data, K, nstart=1, iter.max=10)$cluster
  W   <- sapply(seq(K), function(x) as.numeric(cl==x))
  return(W)
}

# generate methods settings 
B <- mset_gmix(K = c(2,3), erc = c(1,10), init=c(compute_init, "kmed"))


# select setup 2: K=2, erc=10, init = compute_init
mb2  <- B[[2]]
fit2 <- mb2$fn(dat)   
fit2
}

