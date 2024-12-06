\name{plot.bqs}
\alias{plot.bqs}
\title{Plot (Bootstrap) Quadratic Score Results}
\description{
  Produce a plot of bqs (Bootstrap Quadratic Scores). This function creates plots based on the BQS (Bootstrap Quality Scores) data.
}
\usage{
\method{plot}{bqs}(x, score = NULL, perc_scale = FALSE, top = NULL, annotate = NULL, ...)
}
\arguments{
  \item{x}{
    An S3 object of class \code{bqs} as returned by the \code{bqs} function. \code{x} is expected to have the component \code{rankby} set.
  }
  \item{score}{
    Character vector specifying the score(s) to be plotted. Valid scores are \code{"hard"}, \code{"smooth"}, \code{"oob_hard"}, and \code{"oob_smooth"}. If \code{NULL} (default), all valid scores present in \code{x} are plotted.
  }
  \item{perc_scale}{
    Logical; if \code{TRUE}, scales the plot using percentages, relative to the best score. Default is \code{FALSE}.
  }
  \item{top}{
    Numeric; specifies the number of top models to individually highlight. Must be a single number less than or equal to the length of \code{x$methodset}. If \code{NULL} (default), \code{top} is automatically determined based on the \code{score} values.
  }
  \item{annotate}{
    Logical; if \code{TRUE}, annotates the top models in the plot. Default is automatically determined (\code{TRUE} if the number of methods \code{M <= 30}, \code{FALSE} otherwise).
  }
  \item{...}{
    Further arguments passed to or from other methods.
  }
}
\value{
  A plot displaying the Bootstrap Quality Scores.
}
\seealso{
  \code{\link{bqs}}
}
\examples{
\donttest{
# load data
data("banknote")
dat <- banknote[-1]

# set up methods
mlist <- mset_gmix(K=1:3, erc=c(1,100))

# perform bootstrap
# change B and ncores to a much larger value in real problems
res <- bqs(dat, mlist, B = 3, type = "both", rankby="lq",
         ncores = 1, oob = TRUE, savescores = FALSE, saveparams = FALSE)

# Plot with default settings
plot(res)

# Plot in percentage scale relative to first model
plot(res, perc_scale = TRUE)
}
}
