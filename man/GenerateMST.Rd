% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{GenerateMST}
\alias{GenerateMST}
\title{Generate MSTs from clustering results}
\usage{
GenerateMST(X, clusterings)
}
\arguments{
\item{X}{a matrix that specifies the embedding (dim.reduced data).}

\item{clusterings}{a named list of clustering results.}
}
\value{
a list of MSTs
}
\description{
The function generates MSTs using the Slingshot distance method.
}
\keyword{MST}
\keyword{ape}
