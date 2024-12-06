\name{bqs}
\alias{bqs}
\alias{qcluster}

\title{Bootstrapping quadratic scores}

\description{Estimates the expected quadratic score for clustering
  solutions provided by a list of of candidate model, methods or
  algorithmic setting.}


\usage{
   bqs(data,
       methodset,
       B = 10,
       type = "smooth",
       oob = FALSE, 
       ncores = detectCores() - 2,
       alpha = 0.05,
       rankby = ifelse(B == 0, "mean", "lq"),
       boot_na_share = 0.25,
       savescores = FALSE,
       saveparams = FALSE)
}

\arguments{
  \item{data}{
    a numeric vector, matrix, or data frame of observations. Rows
    correspond to observations and columns correspond to
    variables/features.
    Categorical variables and \code{NA} values are not allowed.
  }
  \item{methodset}{
    a list of functions. A function in the list  takes \code{data} as
    input and provide a clustering method to be scored
    (see \emph{Details}).
  }
  \item{B}{
    a integer \code{>=0}.
    If \code{B=0}, the funciton fits and scores the clustering methods on
    the entire data set without resampling (see \emph{Details}).
    \code{B>=1}, sets the number of boostrap replicates (see
    \emph{Details}).
  }
  \item{type}{
      character string specifying the type of score,
      Possibile values are \code{{"smooth", "hard", "both"}}. 
      If \code{="smooth"} (default), only the smooth score is
      estimated.
      If \code{="hard"}, only the hard score is
      estimated. 
      \code{="both"},  both the smooth and the hard
      scores are estimated.
  }
  \item{oob}{
    logical or character string specifying if out-of-bag bootstrap
    is performed.
    Possibile values are \code{{FALSE, TRUE, "only"}}
    If \code{=FALSE} (default),  out-of-bag boostrap is not performed. 
    If \code{=TRUE}, out-of-bag bootstrap is performed along with the
    empirical bootstrap sampling.
    If \code{="only"}, only the out-of-bag bootstrap is performed. 
  }
  \item{ncores}{
   an integer,   it defines the number of cores used for parallel
   computing (see \emph{Details}).
 }
 \item{alpha}{
   a number in \code{(0,1)}, the confidence-level for empirical
   bootstrap quantiles (both-tails). 
  }
  \item{rankby}{
    character string specifying how the scored solution are ranked.
    Possible values are \code{{"lq", "mean", "1se"}}.
    With \code{="lq"} (default), the solutions are ranked by maximizing the estimated
    lower limit of the of the \code{1-alpha} bootstrap confidence
    intervalfor  the expected score.
    With \code{="mean"}, the solutions are ranked by maximizing the estimated
    expected score. 
    With \code{="1se"}, the solutions are ranked by maximizing the estimated
    lower limit of the confidence interval for the expected score
    whose semi-length is equal to a \emph{standard error}.
    The expected score's \emph{standard error} is approximated using the
    boostrap distribution.
  }
  \item{boot_na_share}{
    a numeric value in \code{(0,1)}.
    During the boostrapping a method's score is set to \code{NA} if the
    underlying comptutation runs into errors. 
    Methods resulting in more than \code{B * boot_na_share} errors are
    excluded from the comparison
  }
  \item{savescores}{
    logical, if \code{=TRUE} it returns estimated scores for each
    boostrap sample. 
  }
  \item{saveparams}{
    logical, if \code{=TRUE} it returns estimated cluster parameters for
    each boostrap sample.
  }
}


\value{
An S3 object of class \code{bqs}. Output components are as follows:


\item{smooth}{
  data.frame returned if
  \code{type="smooth"} or \code{type="both"}.
  It contains a summary of the estimated score. The rows corresponds to
  competing methods in \code{methodset} sorted by the specified ranking
  criterion. 
  Columns are as follows:
  \describe{
    \item{\code{id:}}{
      index of the method in the corresponding \code{methodset} list;
    } 
    \item{\code{rank:}}{
      rank of the clustering solution according to \code{rankby};
    } 
    \item{\code{mean:}}{
      expected score;}
    \item{ \code{sterr}}{
      standard error for the mean score;
    }
    \item{\code{lower_qnt:}}{
      lower limit of the confidence interval for the mean score;
    }
    \item{\code{upper_qnt:}}{
      upper limit of the confidence interval for the mean score;
    }
    \item{\code{n_obs:}}{
      number of valid bootstrap samples after filtering erroneous cases;
    }
    \item{\code{n_missing:}}{
      number of filtered erroneous cases.
    }
  }
}
\item{hard}{
  data.frame returned if 
  \code{type="hard"} or \code{type="both"}.
  It reports the results about the hard score in analogy to the previous object
  \code{smooth}.
}
\item{obb_smooth}{
  data.frame returned if 
  \code{type="smooth"} or \code{type="both"} and \code{obb=TRUE} or \code{obb="only"}.
  It reports the results about the smooth score estimated using
  out-of-bang-boostrap samples  in analogy to the previous objects
  \code{smooth} and \code{hard}.
}
\item{obb_hard}{
  data.frame returned if
  \code{type="hard"} or \code{type="both"} and   \code{obb=TRUE} or \code{obb="only"}.
  It reports the results about the hard score estimated using
  out-of-bang-boostrap samples  in analogy to the previous objects
  \code{smooth} and \code{hard}.
}
%%
\item{best_smooth}{
  Clustering produced by the \emph{best} method, according the
  specified \code{rankby} criterion, applied to the smooth score
  estimated using the empirical bootstrap sampling.
}
\item{best_hard}{  
  clustering produced by the \emph{best} method, according the
  specified \code{rankby} criterion, applied to the hard score
  estimated using the empirical bootstrap sampling.
}
\item{best_obb_smooth}{
  clustering produced by the \emph{best} method, according the
  specified \code{rankby} criterion, applied to the smooth score
  estimated using the out-of-bag-boostrap samples.
}
\item{best_obb_hard}{
  clustering produced by the \emph{best} method, according the
  specified \code{rankby} criterion, applied to the hard score
  estimated using the out-of-bag-boostrap samples.
}
%%
\item{data}{
  a list containing information about the input \code{data} set
  necessary for the fruition of the returned object.
}
\item{B}{
  number of bootstrap replicates.
}
\item{methodset}{
  The elements of \code{methodset} for which a solution is
  produced.
}
\item{rankby}{
  the ranking criterion.
}
\item{raw}{
  a list that  allows tracing the bootstrap sampling in almost
  every stage.
  Let \code{n=}sample size, \code{B=}bootstrap samples, \code{M=} number
  of methods in \code{methoset}.
  It contains the following objects.
  \describe{
    \item{\code{boot_id}:}{
      an array of dimension \code{n x B} where
      the \code{j}-th column contains the indexes of the observed data
      points belonging to the \code{j}-th bootstrap sample. That is,  
      \code{data[boot_id[ ,j], ]} gives the \code{j}-th bootstrap data
      set.}
    %%
    \item{\code{scores}:}{
      an array of dimension \code{(M x 3 x B)} returned if
      \code{savescores=TRUE}.
      It reports both hard and smooth scores estimated in each  bootstrap
      replicate.
      \code{score[,1,]} reports a \code{code=1} if the corresponding
      bootstrap sample has been excluded because of errors (otherwise
      \code{code=0}).}
    %% 
    \item{\code{oob_scores}:}{
      returned if  \code{obb=TRUE} or \code{obb="only"} and \code{savescores=TRUE}.
      It is an array is organized as the previous object \code{score}
      but contains information about out-of-bag-bootstrap estimates.
    }
    %%
    \item{\code{params}:}{
      a list returned if \code{saveparams=TRUE}.  \code{params[[m]]}
      contains  estimated cluster parameters for \code{methodset[[m]]}
      where \code{m=1,...,M}.  Each member of the list is a list of
      length \code{B}  where \code{params[[m]][[b]]} contains the
      cluster  parameters fitted by the \code{m}-th method on the
      \code{b}-th  bootstrap sample.
   }
  }
 }
}



\details{
   The function implements the estimation and selection of an
   appropriate  clustering based on the methodology proposed in Coraggio
   and Coretto (2023). In addition, we add the possibility of obtaining
   score estimates using out-of-bag-boorstrap sampling alongside the
   empirical  bootstrap-based estimates proposed in the aforementioned
   paper. Note that the out-of-bag-boorstrap estimates are obtained
   using the same samples used for the emprical bootsrap, therefore,
   \code{oob=TRUE} add a small computational cost.

   
   \strong{Choice of \code{B}.} In theory \code{B} should be as large as
   possible, however, if the list of methods is large and the
   computational capacity is modest, a large \code{B} may require 
   long run times. Coraggio and Coretto (2023) show experiments where
   changing from \code{B=1000} to \code{B=100} introduces a marginal
   increase  in variability. \code{B=100} should be considered as a
   lower  bound. In the case where one has very large method lists,
   high-dimensional  datasets and demanding methods, a possible strategy
   to reduce the computational cost is as follows:
   \enumerate{
     \item set a small value of \code{B}, e.g., \code{B=50} or even less.
     \item Analyze the methods' ranking and identify those methods that
     report score values that are small compared to the top performers.
     \item Narrow down the \code{methodset} list and repeat the
     bootstrap  estimation with a value of \code{B} that is as large as
     possible relative to available computational resources.
   }

   \strong{Parallel computing.} Bootstrap sampling is performed using
   \code{foreach}-based parallel computation via the \code{doParallel}
   parallel backend.   Note that depending on the system settings, and
   how the functions in \code{methodset} make use of parallelism and/or
   multi-threading computing, increasing \code{ncores} may not produce
   the desired reduction in computing time. For instance, this happens
   when  using linear algebra routines optimized for multi-threaded
   computing  (e.g., OpenBLAS, Intel Math Kernel Library (MKL), and so
   on).  These optimized shared libraries already implement
   multi-threading,  and it is necessary to find the optimal trade-off
   between distributing processes over physical cores and
   multi-threading  at the logical unit level of the same physical unit.
   Unfortunately, there is no universal recipe and much depends on
   hardware  architectures, operating system, shared libraries, etc.
   We obtained the best results using OpenBLAS (the tagged
   \emph{serial}) and setting \code{ncores=} the number of physical
   cores. 

   \strong{\code{methodset} argument.} The \code{methodset} argument 
   allows in input a \code{function}, \code{list}, or output from mset 
   functions: \code{mset_user}, \code{mset_gmix}, \code{mset_kmeans}, 
   \code{mset_pam}. It is also possible give any combination of these, 
   concatenated with the \code{mbind} function. When passing a 
   \code{function}, either as a single element or in a list, this must 
   take the data set as its first argument, and must return in output 
   at least \code{list} named \emph{"params"}, conforming with the 
   return value of \code{clust2params}, i.e. a list containing 
   \code{proportion}, \code{mean} and \code{cov} elements, representing 
   the estimated clusters' parameters.
}


\references{
  Coraggio, Luca and Pietro Coretto (2023). Selecting the number of
  clusters, clustering models, and algorithms. A unifying approach based
  on the  quadratic discriminant score.
  \emph{Journal of Multivariate Analysis}, Vol. 196(105181), 1-20.
  doi: \doi{10.1016/j.jmva.2023.105181}
}



\seealso{
  \code{\link{mset_user}}, \code{\link{mset_gmix}}, \code{\link{mset_kmeans}}, 
  \code{\link[cluster]{pam}}, \code{\link{mbind}}, \code{\link{clust2params}}
}


\examples{

# load data
data("banknote")
dat <- banknote[-1]

## set up methods
## see also help('mset_user')  and related functions 
KM   <- mset_kmeans(K = 3)
GMIX <- mset_gmix(K=3, erc=c(1,100)) 

# combine tuned methods
mlist <- mbind(KM, GMIX)

# perform bootstrap
# change B and ncores to a much larger value in real problems 
res <- bqs(dat, mlist, B = 3, type = "both", rankby="lq",
         ncores = 1, oob = TRUE, savescores = TRUE, saveparams = FALSE)
res


\dontrun{
# The following example is more realistic but may take time
# ----------------------------------------------------------
# load data
data("banknote")
dat <- banknote[-1]

# set up kmeans, see help('mset_kmeans')
KM    <- mset_kmeans(K = 2:5)

# set up Gaussian model-based clustering via gmix() 
GMIX  <- mset_gmix(K=2:5, erc=c(1, 50 , 100)) 

# set up Gaussian model-based clustering via library("mclust")
# see examples in help('mset_user')
require(mclust)
mc_wrapper <- function(data, K, ...){
    y <- Mclust(data, G = K, ...)
    y[["params"]] <- list(proportion = y$parameters$pro,
                          mean = y$parameters$mean,
                          cov = y$parameters$variance$sigma)
    return(y)
    }
MC <- mset_user(fname = "mc_wrapper", K = 2:5, modelNames = c("EEI", "VVV"))


# combine tuned methods
mlist <- mbind(KM, GMIX, MC)

# perform bootstrap
# set 'ncores' to the number of available physical cores
res <- bqs(dat, mlist, B = 100, type = "both", rankby="lq", ncores=1,
           oob = TRUE, savescores = TRUE, saveparams = FALSE)
res
}  
}


