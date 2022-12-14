% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\name{VizCellConnectivity}
\alias{VizCellConnectivity}
\alias{VizCellConnectivity.SingleCellExperiment}
\alias{VizCellConnectivity,SingleCellExperiment-method}
\title{Visualize the cell connectivity measure}
\usage{
VizCellConnectivity.SingleCellExperiment(object, dim.red.type, viz.dim.red)

\S4method{VizCellConnectivity}{SingleCellExperiment}(object, dim.red.type = NULL, viz.dim.red = NULL)
}
\arguments{
\item{object}{an object of class \code{SingleCellExperiment}.}

\item{dim.red.type}{a character that specifies the name of the dimensionality
reduction method saved in the \code{SingleCellExperiment} object.}

\item{viz.dim.red}{a matrix for the user-provided, custom dimensionality
reduction.}
}
\value{
a \code{ggplot2} object
}
\description{
The function visualizes the cell connectivity measure, which is used as 
a basis for determining the shape of the trajectory. A low value indicates 
that the cell population has, on average, a small number of connections 
with other cell types. These are generally the leaf/end nodes of 
the trajectory tree. Cell subsets that have a high value are 
at the branching points of the trajectory, and they have, on average, more
connections compared to the other cell types. The measure is calculated by
counting in each MST how many edges the clusters have, divided by the number
of clusters in each clustering, and by averaging the ratios of all the MSTs.
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
dim_red <- dyndimred::dimred_mds(t(logcounts(sce)),ndim=2)
VizCellConnectivity(sce,viz.dim.red=dim_red)
}
\keyword{cell}
\keyword{connectivity}
\keyword{visualization}
