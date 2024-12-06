\name{qscore}

\alias{qscore}

\title{Clustering Quadratic Score}

\description{
   Computes both the hard and the smooth quadratic score of a clustering.
   Handles both 
}

\usage{
   qscore(data, params, type = "both")
   }

\arguments{
   \item{data}{a numeric vector, matrix, or data frame of observations. 
   Rows correspond to observations and columns correspond to variables/features. 
   Let \code{N=nrows(data)} and \code{P=ncol(data)}. 
   Categorical variables and \code{NA} values are not allowed.
 }

  \item{params}{
   a list containing  cluster parameters \emph{(size, mean, cov)}. Let 
   \code{K=}\emph{number of clusters}.
   The elements  of the list are as follows: 
   \code{$prop=}vector of clusters' proportions; 
   \code{$mean=}matrix of dimension \code{(P x K)} containing the clusters' mean 
   parameters;
   \code{$cov=}array of size \code{(P x P x K)} containing the clusters'
   covariance matrices.
   }
   
  \item{type}{the type of score, a character in the set 
    \code{c("both", "smooth", "hard")}. The default value is set to
    "booth". See \emph{Details}.
  }
}


\value{
A numeric vector with both the hard and the smooth score, or only one of them 
depending on the argument \code{type}.
}


\details{
The function calculates quadratic scores as defined in equation (22) in Coraggio
and Coretto (2023). 
}


\references{
  Coraggio, Luca and Pietro Coretto (2023). Selecting the number of
  clusters, clustering models, and algorithms. A unifying approach based
  on the  quadratic discriminant score.
  \emph{Journal of Multivariate Analysis}, Vol. 196(105181), 1-20.
  DOI: \doi{10.1016/j.jmva.2023.105181}
  }

\seealso{
  \code{\link{clust2params}}
}


\examples{
# --- load and split data
data("banknote")
set.seed(345)
idx   <- sample(1:nrow(banknote), size = 25, replace = FALSE)
dat_f <- banknote[-idx, -1] ## training data set  
dat_v <- banknote[ idx, -1] ## validation data set 


# --- Gaussian model-based clustering, K=3
# fit clusters 
fit1 <- gmix(dat_f, K=3)
## compute quadratic scores using fitted mixture parameters
s1 <- qscore(dat_v , params = fit1$params)
s1


# --- k-means clustering, K=3
# obtain the k-means partition 
cl_km <- kmeans(dat_f, centers = 3, nstart = 1)$cluster
## convert k-means hard assignment into  cluster parameters 
par_km <- clust2params(dat_f, cl_km)
# compute quadratic scores 
s2 <- qscore(dat_v, params = par_km)
s2
}


