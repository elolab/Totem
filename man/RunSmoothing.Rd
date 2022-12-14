% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\name{RunSmoothing}
\alias{RunSmoothing}
\alias{RunSmoothing.SingleCellExperiment}
\alias{RunSmoothing,SingleCellExperiment-method}
\title{Run Slingshot smoothing}
\usage{
RunSmoothing.SingleCellExperiment(object, slingshot.par.list)

\S4method{RunSmoothing}{SingleCellExperiment}(
  object,
  slingshot.par.list = list(shrink = 1, reweight = TRUE, reassign = TRUE, thresh = 0.001,
    maxit = 10, stretch = 2, smoother = "smooth.spline", shrink.method = "cosine",
    dist.method = "slingshot")
)
}
\arguments{
\item{object}{an object of class \code{SingleCellExperiment}.}

\item{slingshot.par.list}{a list of the Slingshot parameters.}
}
\value{
an object of class \code{SingleCellExperiment}
}
\description{
The function runs the Slingshot smoothing. The Slingshot parameters can also
be adjusted. The Slignshot smoothing is performed for the selected 
clusterings that were selected using \code{SelectClusterings} function.
}
\examples{
## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
library(SingleCellExperiment)
data(binary_tree_1)
sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
logcounts = t(binary_tree_1$normalized_data)))
sce <- PrepareTotem(sce)
sce <- RunDimRed(sce,dim.red.method="lmds",
dim.reduction.par.list=list(ndim=5))
sce <- RunClustering(sce,N.clusterings=200) #Use default N.clusterings=10000
sce <- SelectClusterings(sce)
sce <- RunSmoothing(sce)
}
\keyword{slingshot}
\keyword{smoothing}
