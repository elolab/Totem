% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\name{ReturnSlingshotObject}
\alias{ReturnSlingshotObject}
\alias{ReturnSlingshotObject.SingleCellExperiment}
\alias{ReturnSlingshotObject,SingleCellExperiment-method}
\title{Return a trajectory as a PseudotimeOrdering object}
\usage{
ReturnSlingshotObject.SingleCellExperiment(object, traj.name)

\S4method{ReturnSlingshotObject}{SingleCellExperiment}(object, traj.name)
}
\arguments{
\item{object}{an object of class \code{SingleCellExperiment}}

\item{traj.name}{a character for the trajectory name}
}
\value{
an object of class \code{PseudotimeOrdering}
}
\description{
The function returns a smoothed trajectory as a \code{PseudotimeOrdering} 
object, which is the object class of Slingshot.
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
sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
sce <- SelectClusterings(sce,selection.method = 1,selection.N.models = 2)
sce <- RunSmoothing(sce)
clustering_name <-  ReturnTrajNames(sce)[1]
dynwrap_object <- ReturnSlingshotObject(sce,clustering_name)
}
\keyword{PseudotimeOrdering}
\keyword{Slingshot}
\keyword{object}
\keyword{smooth}
\keyword{trajectory}
