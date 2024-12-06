\name{predict.mbcfit}
\alias{predict.mbcfit}

\title{Predict Hard Clustering Assignments for using Mixture Models}

\description{
  This function predicts cluster assignments for new data based on an existing model of class \code{mbcfit}. The prediction leverages information from the fitted model to categorize new observations into clusters.
}

\usage{
\method{predict}{mbcfit}(object, newdata, ...)
}

\arguments{
  \item{object}{
    An object of class \code{mbcfit}, representing the fitted mixture model. This is typically obtained in output from the \code{gmix} function. See \emph{Details}.
  }
  \item{newdata}{
    A numeric vector, matrix, or data frame of observations. Rows correspond to observations and columns correspond to variables/features. Categorical variables and \code{NA} values are not allowed. The number of columns must be coherent with that implied by \code{x}. See \emph{Details}.
  }
\item{...}{
  Further arguments passed to or from other methods.
}
}

\value{
  A vector of length \code{nrow(data)} containing the estimated cluster labels for each observation in the provided \code{data}.
}

\details{
  The \code{predict.mbcfit} function utilizes the parameters of a previously fitted \code{mbcfit} model to allocate new data points to estimated clusters. The function performs necessary checks to ensure the \code{mbcfit} model returns valid estimates and the dimensionality of the new data aligns with the model.

  The \code{mbcfit} object must contain a component named \code{params}, which is itself a list containing the following necessary elements, for a mixture model with K components:
    \describe{
        \item{\code{proportions}}{A numeric vector of length K, with elements summing to 1, representing cluster proportions.}
        \item{\code{mean}}{A numeric matrix of dimensions \code{c(P, K)}, representing cluster centers.}
        \item{\code{cov}}{A numeric array of dimensions \code{c(P, P, K)}, representing cluster covariance matrices.}
    }

  Data dimensionality is \code{P}, and new data dimensionality must match (\code{ncol(data)} must be equal to \code{P}) or otherwise the function terminates with an error message.

  The predicted clustering is obtained as the MAP estimator using posterior weights of a Gaussian mixture model parametrized at \code{params}.
  Denoting with \eqn{z(x)}{z(x)} the predicted cluster label for point \eqn{x}{x}, and with \eqn{\phi}{phi} the (multivariate) Gaussian density:
  \deqn{z(x) = \underset{k=\{1,\ldots,K\}}{\arg\,\max} \frac{\pi_k\phi(x, \mu_k, \Sigma_k)}{\sum_{j=1}^K \pi_j\phi(x, \mu_j, \Sigma_j)}}{z(x) = argmax_k { (pi_k * phi(x, mu_k, Sigma_k)) / (sum_{j=1}^K (pi_j * phi(x, mu_j, Sigma_j))) }}
}

\references{
  Coraggio, Luca and Pietro Coretto (2023). Selecting the number of
  clusters, clustering models, and algorithms. A unifying approach based
  on the  quadratic discriminant score.
  \emph{Journal of Multivariate Analysis}, Vol. 196(105181), 1-20.
  doi: \doi{10.1016/j.jmva.2023.105181}
}

\seealso{
  \code{\link{gmix}}
}

\examples{
# load data
data(banknote)
dat <- banknote[,-1]

# Estimate 3-components gaussian mixture model
set.seed(123)
res <- gmix(dat, K = 3)

# Cluster in output from gmix
print(res$cluster)

# Predict cluster on a single point
# (keep table dimension)
predict(res, dat[1, , drop=FALSE])

# Predict cluster on a subset
predict(res, dat[1:10, ])

# Predicted cluster on original dataset are equal to the clustering from the gmix model
all(predict(res, dat) == res$cluster)
}