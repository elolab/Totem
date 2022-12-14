% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\name{SelectClusterings}
\alias{SelectClusterings}
\alias{SelectClusterings.SingleCellExperiment}
\alias{SelectClusterings,SingleCellExperiment-method}
\title{Perform clustering selection}
\usage{
SelectClusterings.SingleCellExperiment(
  object,
  selection.method,
  selection.N.models,
  selection.k.range,
  selection.stratified,
  prior.clustering
)

\S4method{SelectClusterings}{SingleCellExperiment}(
  object,
  selection.method = 1,
  selection.N.models = 1,
  selection.k.range = NULL,
  selection.stratified = FALSE,
  prior.clustering = NULL
)
}
\arguments{
\item{object}{an object of class \code{SingleCellExperiment}}

\item{selection.method}{a positive integer specifying 
the clustering selection method.
\code{1} will use the Variance Ratio Criterion (Calinski-Harabasz score) 
VRC with the cell connectivity measure as input.
\code{2} will use the VRC with the dimensionally reduced data as input.
\code{3} will use the mean of \code{1} and \code{2}.
\code{4} will use the mean of \code{1} and the ARI between a clustering and 
the \code{prior.clustering}. ARI is the adjusted Rand index that calculates 
the similarity of two clusterings. \code{prior.clustering} parameter 
specifies the prior clustering that the user wants to use to direct the 
clustering search towards clusterings that are more similar with 
\code{prior.clustering}. \code{5} will use the mean of \code{2} and 
the ARI between a clustering and 
the \code{prior.clustering}. ARI is the adjusted Rand index that measures the
similarity of two clusterings.\code{prior.clustering} parameter specifies the
prior clustering that the user wants to use to direct the clustering search
towards clusterings that are more similar with \code{prior.clustering}.}

\item{selection.N.models}{a positive integer specifying the number of
top-ranking clusterings that are selected. If the value is a a numeric 
vector, this selects a set of clusterings of which position in the ranked
list matches the integers (e.g., \code{2:4} would select 
the clusterings ranked 2nd,3rd and 4th).}

\item{selection.k.range}{an integer vector that specifies 
the range of the number of clusters that are allowed in the selection. E.g.
\code{c(10,11)} would select clusterings that have 10 or 11 clusters.}

\item{selection.stratified}{a logical variable (TRUE or FALSE) that can be
used to activate the stratified selection method. If activated, 
the number of clusters will vary evenly in the selected clusterings. Can
be useful if the top results include trajectories that are not satisfactory,
and the user wants to try diverse trajectories.}

\item{prior.clustering}{a vector that includes a prior clustering. When used
in conjunction with \code{selection.method=4} or \code{selection.method=5},
it can be used to select clusterings that are more similar with
the prior clustering.}
}
\value{
an object of class \code{SingleCellExperiment}
}
\description{
The function selects clusterings that are used in the 
trajectory smoothing with Slingshot. The user can choose from
different criteria to perform the selection.
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
sce <- SelectClusterings(sce)
}
\keyword{ARI}
\keyword{clustering}
\keyword{selection}
