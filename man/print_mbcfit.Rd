\name{print.mbcfit}
\alias{print.mbcfit}

\title{Display Information for Mixture Model Objects}

\description{
  This function provides a print method for objects of class \code{mbcfit}, returned in output by the \code{gmix} function.
}

\usage{
\method{print}{mbcfit}(x, ...)
}

\arguments{
  \item{x}{
    An object of class \code{mbcfit}, typically a result of the \code{gmix} function.
  }
  \item{...}{
    Further arguments passed to or from other methods.
  }
}
\value{No return value, called for side effects}
\details{
  The \code{print.mbcfit} function gives a summary of a model-based clustering fit, estimated using the \code{gmix} function.

  The function handles different \code{code} values from the object's \code{info} field, each representing a specific status or error condition:
  \describe{
    \item{\code{-2}}{
      'Lapack DSYEV failed'. This error occurs whenever any of the cluster-covariance matrices becomes singular during estimation, using the EM algorithm.
    }
    \item{\code{-1}}{
      'Memory allocation error'. This error occurs when there is insufficient available memory to allocate the quantities required to execute the EM algorithm.
    }
    \item{\code{1}}{
      Success.
    }
    \item{\code{2}}{
      'gmix' did not converge (iterations reached the maximum limit).
    }
    \item{\code{3}}{
      EM algorithm failed; no better than the initial solution. This error occurs whenever the EM algorithm failed for other reasons (e.g., degenerate posterior-weights could not be prevented), and it was not possible to find a solution.
    }
  }

  The printed output also lists available components of the \code{mbcfit} object and summarizes the number of clusters found and their size, whenever this information is available.
}

\seealso{
  \code{\link{gmix}}
}

\examples{
set.seed(123)

# Estimate a simple a 3-clusters Gaussian mixture model, using iris data as example
res <- gmix(iris[,-5], K = 3, erc = 10)

# Print the 'gmix' output
print(res)
}