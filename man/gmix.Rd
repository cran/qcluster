\name{gmix}

\alias{gmix}

\title{Gaussian Mixture Modelling}

\description{
Fast implementation of the EM algorithm for ML estimation and clustering of Gaussian mixture models with covariance matrix regularization based on eigenvalue 
ratio constraints.
}


\usage{
gmix(
   data,
   K = NA,
   erc = 50,
   iter.max = 1000,
   tol = 1e-8,
   init = "kmed",
   init.nstart = 25,
   init.iter.max = 30,
   init.tol = tol,
   save_cluster = TRUE,
   save_params = TRUE,
   save_taus = FALSE)
}

\arguments{

  \item{data}{
    a numeric vector, matrix, or data frame of observations. Rows
    correspond to observations and columns correspond to
    variables/features.  Let \code{N=nrows(data)} and
    \code{P=ncol(data)}. Categorical variables and \code{NA} values are
    not allowed.
  }

  \item{K}{
    the number of mixture components or clusters. It can be left
    \code{NA} for certain specifications of \code{init} (see below) where
    the number of cluster is retrieved from the initial partition.
  }

  \item{erc}{
    a numeric value \code{>=1} specifying the eigenvalue ratio
    constraint (See \emph{Details}).
  }

  \item{iter.max}{
    maximum number of iterations for the EM algorithm.
  }

  \item{tol}{
    tolerance for the convergence of the EM algorithm.
  }

  \item{init}{
    a character in the set \code{c("kmed", "kmeans", "pam")}, a vector, a matrix, 
    or a callable giving the initial assignment of data points (see
    \emph{Details}). The default choice is "kmed".
  }

  \item{init.nstart}{
    number of initial partitions ((see \emph{Details})).
  }

  \item{init.iter.max}{
    maximum number of iterations for each run of the kmedian
    initialization.
  }

  \item{init.tol}{
    tolerance for the convergence of each ran of the  kmedian
    initialization.
  }

  \item{save_cluster}{
    logical, if \code{TRUE} the point-to-cluster assignment based on the 
    \emph{maximum a posteriori probability} (MAP) rule is returned.
  }

  \item{save_params}{
    logical, if \code{TRUE} the estimated mixture parameters are
    returned.
  }

  \item{save_taus}{
    logical, if \code{TRUE} the posterior class probabilities are returned
    (these are also known as \emph{posterior weights} or \emph{fuzzy
      weights}).
  }
}





\value{
An S3 object of class \code{'mbcfit'}. Output components are as follows:

\item{info}{
   information on convergence and errors (see notes).
   }
\item{iter}{
   number of iterations performed in the underlying  EM-algorithm.
   }
\item{N}{
   number of data points.
}
\item{P}{
   data dimension.
   }
\item{K}{
   number of clusters.
   }
\item{eloglik}{
   sample expected log\code{-}likelihood.
   }
\item{size}{
   cluster size (counts).
   }
\item{cluster}{
   cluster assignment based on the  \emph{maximum a posteriori} rule  (MAP).
   }
\item{taus}{
   a matrix of dimension \code{(N x {K})} where \code{tau[i, k]} is the 
   estimated posterior probability that the \emph{i}th observation belongs to the
   \emph{k}th cluster.
   }
\item{params}{
   a list containing mixture components parameters.   
   The elements  of the list are as follows:  
   \code{$prop=}vector of proportions;  
   \code{$mean=}=matrix of dimension \code{(P x K)} containing mean parameters; 
   \code{$cov=}array of size \code{(P x P x K)} containing covariance 
   matrices. 
 }
 \item{info}{
   a list with two components named giving information about underlysing
   EM algorithm. The  \code{code} obejects can take the following values:
   \itemize{
     \item \code{code=1}: the algorithm converged within \code{iter.max}. 
     \item \code{code=2}: the algorithm reached \code{iter.max}.
     \item \code{code=3}: the algorithm did not move from initial values.
     \item \code{code=-1}: unexpected memory allocation issues occured.
     \item \code{code=-2}: unexpected LAPACK routines errors  occured.
   }
   The  \code{flag} obejects can take the following values:
   \itemize{
     \item \code{flag=0} no flag.
     \item \code{flag=1} numerically degenerate posterior probabilities.
     (\code{taus}) could not be prevented.
     \item \code{flag=2} the ERC was enforced at least once.
     \item  \code{flag=3} if condition of \code{flag=1} and \code{flag=2}
     occurred.
   }
 }
}



\details{
  The function implements the constrained ML estimator studied in
  Coretto and Hennig  (2023). The convariance matrix constraints are
  computed according to the \code{CM1-step} of  Algorithm 2 of Coretto
  and Hennig (2017). This function uses highly optimized C code for fast
  execution. The constrained M-step extensively uses low-level common
  linear algebra matrix operations (BLAS/LAPACK routines). Consequently,
  to maximize computational efficiency, it is recommended that the best
  available shared libraries, such as  OpenBLAS, Intel Math Kernel
  Library (MKL),  etc., be set up. 

  \strong{Initialization.} 
  The default method, set with \code{init="kmed"}, uses fast C
  implementation  of the  k\code{-}medians algorithm with random initial
  centers drawn uniformly over the \code{data} rows \code{init.iter.max}
  times. Depending on the computer power available
  it is suggested to set \code{init.iter.max} as large as  possible
  particularly  in cases where the data set dimensionality is large in
  terms  of both sample  size and number of  features.
  Setting  \code{init="kmeans"} one replaces the K\code{-}medians with
  the  K\code{-}means. With \code{init="pam"}  initial clusters are
  determined  using the PAM algorithm based on  Euclidian distances. The
  latter  does not perform multiple starts. 
  The user can also set \code{init = x} where \code{x} is a vector of
  integers  of length \code{N=nrow(data)} representing an initial hard
  assignment  of data points to the mixture components or clusters (see \emph{Examples}). 
  Another possibility is to set \code{init = W} where \code{W} is a
  matrix of  dimension \code{(N x {K})} containing the initial posterior
  probabilities that the \emph{i}th observation belongs to the
  \emph{k}th  cluster. The assignment provided via \code{W} can be hard
  (0\code{-}1 weights with the constraint that only a 1 is possible in
  each row of \code{W}, or smooth (each row of \code{W} must sum up to
  1).   \code{W} can be seen as the initial version of the object
  \code{tau}  describied in the \emph{Value} section below.
  The last alternative is to set \code{init = f(data)}. Here
  \code{f(data)}  is a function 
  that takes \code{data} as an input and returns the matrix with an
  initial  hard/smooth assignment as the \code{W} matrix previously
  described (see the example below).	


  \strong{Eigenvalue ratio constraint (\code{erc}).} 
  It is  the maximum allowed ratio between within-cluster covariance
  matrix  eigenvalues. 
  It defines the so\code{-}called \emph{eigenratio constraint}. 
  \code{erc=1} enforces spherical clusters with equal covariance matrices. 
  A large \code{erc} allows for large between-cluster covariance discrepancies. 
  It is suggested to never set \code{erc} arbitrarily large, its main role is to 
  prevent degenerate covariance parameters and the related emergence of spurious 
  clusters (see \emph{Referenceses} below).
  Finally, in order to facilitate the setting of \code{erc}, it is
  suggested to scale the  columns of \code{data} whenever measurement
  units  of the different variables 
  are grossly incompatible. 

}



\references{
  Coretto, Pietro and Christian Hennig (2017).
  Consistency, breakdown robustness, and algorithms for robust improper
  maximum  likelihood clustering.
  \emph{Journal of Machine Learning Research}, Vol. 18(142), pp. 1-39.   URL: 
  \href{https://jmlr.org/papers/v18/16-382.html}{https://jmlr.org/papers/v18/16-382.html}
  
  Coretto, Pietro and Christian Hennig (2023)
  Nonparametric consistency for maximum likelihood estimation and clustering based 
  on mixtures of elliptically-symmetric distributions.
  \emph{arXiv:2311.06108}. URL: 
  \href{https://arxiv.org/abs/2311.06108}{https://arxiv.org/abs/2311.06108}
}



\examples{
# --- load data
data("banknote")
dat <- banknote[-1]
n   <- nrow(dat) #sample size
nc  <- 2         #number of clusters


   
# fit 2 clusters using the default k-median initialization 
# In real applications set 'init.nstart' as large as possibile
set.seed(101)
fit1 <- gmix(dat, K = nc, init.nstart = 1)
print(fit1)

# plot partition (default)
plot(x = fit1, data = dat)


# plot partition onto the first 3 principal component coordinates 
plot(x = fit1, data = prcomp(dat)$x, margins = c(1,2,3),
     pch_cl = c("A", "B"), col_cl = c("#4285F4", "#0F9D58"),
     main = "Principal Components")


 
# user-defined random initialization with hard assignment labels 
set.seed(102)
i2   <- sample(1:nc, size = n, replace = TRUE)
fit2 <- gmix(dat, K = 2, init = i2)
plot(x=fit2, data = dat)



# user-defined smooth "toy" initialization: 
# 50% of the points are assigned to cluster 1 with probability 0.95 and to
# cluster 2 with probability 5%. The remaining data points are assigned to
# cluster 1 with probability 10% and  to cluster 2 with probability 10%
# 
set.seed(103)
idx        <- sample(c(TRUE, FALSE), size = n, replace = TRUE)
i3         <- matrix(0, nrow = n, ncol = nc) 
i3[idx,  ] <- c(0.9, 0.1)
i3[!idx, ] <- c(0.1, 0.9)
# fit
fit3  <- gmix(dat, K = nc, init = i3)
plot(x=fit3, data = dat)



# user-defined function for initialization
# this one produces a 0-1 hard posterior matrix W based on kmeans
#
compute_init <- function(data, K){
  cl  <- kmeans(data, K, nstart=1, iter.max=10)$cluster
  W   <- sapply(seq(K), function(x) as.numeric(cl==x))
  return(W)
} 
fit4 <- gmix(dat, K = nc, init = compute_init)
plot(fit4, data = dat)

}


