\name{clust2params}

\alias{clust2params}

\title{Converts Hard Assignment Into Cluster Parameters}

\description{
Transforms cluster labels into a list of parameters  describing cluster size, 
mean, and dispersion.}


\usage{
   clust2params(data, cluster)
}

\arguments{
   \item{data}{
   a numeric vector, matrix, or data frame of observations. 
   Rows correspond to observations and columns correspond to variables/features. 
   Categorical variables and \code{NA} values are not allowed.
   }
   
   \item{cluster}{
   a vector of integers representing cluster labels.
   }
}


\value{
   A list containing  cluster parameters. 
   Let \code{P=}\emph{number of variable/features} and \code{K=}\emph{number of clusters}. 
   The elements  of the list are as follows: 
   \itemize{
   \item \code{prop:} a vector of clusters' proportions; 
   \item \code{mean:} a matrix of dimension \code{(P x K)} containing the clusters' mean 
   parameters;
   \item \code{cov:} an array of size \code{(P x P x K)} containing the clusters'
   covariance matrices.
   }
}
   

\examples{
# load data
data("banknote")

# compute the k-means partition 
set.seed(2024)
cl <- kmeans(banknote[-1], centers = 2, nstart = 1)$cluster

# convert k-means hard assignment into cluster parameters 
clpars <- clust2params(banknote[-1], cl) 
clpars

}


