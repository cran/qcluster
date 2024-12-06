
\name{mset_pam}

\alias{mset_pam}

\title{Generates Methods Settings for Partitioning Around Medoids (Pam) Clustering}

\description{
  The function generates a software abstraction of a list of clustering
  models implemented through the a set of  tuned methods and algorithms. 
  In particular, it generates a list of \code{\link[cluster]{pam}}-type
  functions each combining tuning parameters and other algorithmic settings.  
  The generated functions are ready to be called on the data set.
}

\usage{
mset_pam(K = seq(10),
         metric = "euclidean",
         medoids = if (is.numeric(nstart)) "random",
         nstart = if (variant == "faster") 1 else NA,
         stand = FALSE,
         do.swap = TRUE,
         variant = "original",
         pamonce = FALSE)
}

\arguments{
  \item{K}{
    a vector/list, specifies the number of clusters.  
  }
  \item{metric}{  
   a vector, contains the settings of the \code{metric} parameter of \code{\link[cluster]{pam}}.  
  }
  \item{medoids}{
    list, contains the settings of the \code{medoids} parameter of \code{\link[cluster]{pam}}.  
  }
  \item{nstart}{
    a vector, contains the settings of the \code{nstart} parameter of \code{\link[cluster]{pam}}.  
  }
  \item{stand}{
    a vector, contains the settings  of the \code{stand} parameter of \code{\link[cluster]{pam}}.  
  }
  \item{do.swap}{
    a vector, contains the settings of the \code{do.swap} parameter of \code{\link[cluster]{pam}}.  
  }
  \item{variant}{
    a list, contains the settings of the \code{variant} parameter of \code{\link[cluster]{pam}}.  
  }
  \item{pamonce}{
    a vector, contains the settings of the \code{pamonce} parameter of \code{\link[cluster]{pam}}.  
  }
}


\details{
  The function produces functions implementing competing clustering methods
  based on the PAM clustering methodology  as implemented in
  \code{\link[cluster]{pam}}.
  This is a specialized version of the more general function
  \code{\link{mset_user}}. 
  In particular, it produces a list of \code{\link[cluster]{pam}} functions each
  corresponding to a specific setup in terms of
  hyper-parameters (\emph{e.g.} the number of clusters) and algorithm's
  control parameters (\emph{e.g.}  initialization).
  See \code{\link[cluster]{pam}} for more detail for a detailed description of
  the role of each argument and their data types. 
}



\value{
  An S3 object of class \code{'qcmethod'}. Each element of the list
  represents a competing method containing the following objects
  
  \item{fullname}{
    a string identifying the setup.
  }
  \item{callargs}{
    a list with  \code{\link[cluster]{pam}} function arguments.
  }
  \item{fn}{
   the function implementing the specified setting. This \code{fn}
    function can be executed on the data set.
    It has two arguments: \code{data} and \code{only_params}.
    \code{data} is a data matrix or data.frame
    \code{only_params} is logical.
    If \code{only_params==FALSE} (default), \code{fn} will return the
    object returned by \code{\link[cluster]{pam}}.
    If \code{only_params==TRUE} (default) \code{fn} will return only cluster 
    parameters (proportions, mean, and cov, see \link{clust2params}.
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
  \code{\link[cluster]{pam}},\code{\link{mset_user}}, \code{\link{bqs}}
}







\examples{
# 'pam' settings combining number of clusters K={2,3}, and dissimilarities {euclidean, manhattan}
A <- mset_pam(K = c(2,3), metric = c("euclidean", "manhattan"))
   
# select setup 1: K=2, metric = "euclidean"
m <- A[[1]]
print(m)

      
# cluster with the method set in 'm'
data("banknote")
dat  <- banknote[-1]
fit1 <- m$fn(dat)   
fit1
class(fit1)


# if only cluster parameters are needed
fit1b <- m$fn(dat, only_params = TRUE)   
fit1b

}

