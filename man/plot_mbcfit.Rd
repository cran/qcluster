\name{plot.mbcfit}
\alias{plot.mbcfit}

\title{Plot Fitted Mixture Models}

\description{
  This function provides a plot method for objects of class \code{mbcfit}, returned as output by the \code{gmix} function. It serves as a wrapper around \code{\link{plot_clustering}}, allowing easy visualization of clustering results, including clustering assignments, contours, and boundaries.
}

\usage{
\method{plot}{mbcfit}(x, data = NULL, subset = NULL,
                      what = c("clustering", "contour"),
                      col_cl = NULL, pch_cl = NULL, ...)
}

\arguments{
  \item{x}{
    An object of class \code{mbcfit}, typically a result of the \code{gmix} function.
  }
  \item{data}{
    \code{NULL} or a data matrix, data frame, or array containing data points to be plotted. See \emph{Details}.
  }
  \item{subset}{
    A numeric vector indexing columns of \code{data} to subset and focus the plot on specific features. Default is \code{NULL}.
  }
  \item{what}{
    Character vector specifying which elements to plot. Options are \code{"clustering"}, \code{"contour"}, and \code{"boundary"}. Default is to plot \code{"clustering"} and \code{"boundary"}. See \emph{Details}.
  }
  \item{col_cl}{
    A vector of colors to use for clusters (one for each cluster). Default is \code{NULL}, which uses a default sequence of colors.
  }
  \item{pch_cl}{
    A vector of plotting symbols (one for each cluster) to use for clusters. Default is \code{NULL}, which uses a default sequence of symbols.
  }
  \item{...}{
    Further arguments passed to or from other methods.
  }
}

\details{
  The \code{plot.mbcfit} function provides a plotting method for objects of the class \code{mbcfit}. It acts as a wrapper around the \code{\link{plot_clustering}} function, allowing users to easily generate various plots to analyze the clustering results. A plot is produced only upon a successful \code{mbcfit} estimate, i.e., when \code{mbcfit} has \code{code} equal to either \code{1} or \code{2}.

  When \code{data} is \code{NULL} (the default), the function plots only contour sets (and optionally clustering boundaries) for the estimated mixture density components, using the \code{params} information from the \code{mbcfit} object. When \code{data} is not \code{NULL}, the function additionally plots data points and their hard clustering labels, which are obtained using \code{mbcfit} to predict the cluster labels (see \code{\link{predict.mbcfit}}).

}

\value{
  A plot displaying the data with clustering information, contours, and/or boundaries, depending on the specified \code{what} argument.
}

\seealso{
  \code{\link{gmix}}, \code{\link{plot_clustering}}, \code{link{predict.mbcfit}}
}

\examples{
# load data
data("banknote")
dat <- banknote[-1]

# fit 2 clusters
set.seed(123)
fit <- gmix(dat, K = 2, init.nstart = 1)
print(fit)

# plot partition (default)
plot(x = fit, data = dat)

\donttest{
# plot partition onto the first 3 coordinates
plot(x = fit, data = dat, subset = c(1:3), pch_cl = c("A", "B"),
     col_cl = c("#4285F4", "#0F9D58"), what = "clustering")

# additionally plot clustering boundary and contour sets
plot(x = fit, data = dat, subset = c(1:3), pch_cl = c("A", "B"),
     col_cl = c("#4285F4", "#0F9D58"), what = c("clustering", "boundary", "contour"))
}
}