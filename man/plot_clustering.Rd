\name{plot_clustering}
\alias{plot_clustering}
\title{Plot Data With Clustering Information}
\description{
  This function plots data and optionally adds clustering information such as clustering assignments, contours, or boundaries.
}
\usage{
plot_clustering(data,  subset = NULL,
                cluster = NULL, params = NULL,
                what = c("clustering", "contour", "boundary"),
                col_cl = NULL, pch_cl = NULL)
}
\arguments{
  \item{data}{
    a numeric vector, matrix, or data frame of observations. Rows correspond to observations and columns correspond to variables/features. Categorical variables and \code{NA} values are not allowed.
  }
  \item{subset}{
    A numeric vector indexing columns of \code{data} to subset and focus the plot on specific features. Default is \code{NULL}.
  }
  \item{cluster}{
    A vector of cluster assignments. If provided, the plot can display clustering information as specified in \code{what}. Must have the same number of observations as \code{data}
  }
  \item{params}{
    A list of clustering parameters, including \code{proportion}, \code{mean}, and \code{cov}. If provided, the plot can display contour and boundary information.
  }
  \item{what}{
    Character vector specifying which elements to plot. Options are \code{"clustering"}, \code{"contour"}, and \code{"boundary"}. Default is to plot \code{"clustering"} whenever \code{cluster} is not \code{NULL}, and \code{"contour"} and \code{"boundary"} whenever \code{params} is not \code{NULL}.
  }
  \item{col_cl}{
    A vector of colors to use for clusters (one for each cluster). Default is \code{NULL}, which uses a default sequence of colors.
  }
  \item{pch_cl}{
    A vector of plotting symbols (one for each cluster) to use for clusters. Default is \code{NULL}, which uses a default sequence of symbols.
  }
}
\value{No return value, called for side effects}
\examples{
# Example data
set.seed(123)
data <- rbind(
            matrix(rnorm(100 * 2), ncol = 2),
            matrix(rnorm(100 * 2) + 2, ncol = 2)
        )
cluster <- c(rep(1, 100), rep(2, 100))
params <- clust2params(data, cluster)

# Plot with clustering information
plot_clustering(data, cluster = cluster, what = "clustering")

# Plot with subset of variables
plot_clustering(data, cluster = cluster, subset = 1, what = c("clustering", "contour"))

# Plot with customized colors and symbols
plot_clustering(data, cluster = cluster, params = params,
                col_cl = c("magenta", "orange"), pch_cl = c("A", "B"))
}
\seealso{
  \code{\link{bqs}}, \code{\link{clust2params}}
}