# Totem

A user-friendly tool for inferring tree-shaped trajectories from single cell data. 

## Description

`Totem` package is a tool designed to facilitate the inference of tree-shaped trajectories from single-cell
data [1]. These include linear, bifurcating, multifurcating, and other tree-shaped trajectories.

In trajectory inference methods like `Slingshot`[2] that require a clustering to construct the Minimum Spanning Tree (`MST`) as the backbone of the trajectory, much of the success of the trajectory inference depends on finding a good clustering that produces a biologically meaningful trajectory with the correct milestone connections. The MST construction is determined
by the cluster distances of which generation is sensitive to the cluster structure, the pre-processing steps (dimensionality reduction), and how the distances are calculated.

`Totem` generates a large number of clustering results with a *k*-medoids algorithm (CLARA) and
constructs an MST for each clustering. `Totem` estimates `the connectivity` of the cells (`cell connectivity`)
from the MSTs, which measures how many connections the cells have on average to
other clusters in the MSTs. The cell connectivity allows users to visually see which parts of the embedding are likely to be end points (leaves) and branch points of the trajectory, and the user can leverage the cell connectivity to optimize the trajectory.
`Totem` uses the `Calinski-Harabasz score` and the cell connectivity to rank the clustering results.

From the ranked clusterings the user can select one or several clusterings and easily
visualize the resulting MSTs. This approach is useful because finding a good clustering automatically that will also generate a sensible MST is often difficult, requiring the user to try different tools and parameters. The user can also provide `a prior clustering` to direct the search towards clusterings that are more similar with the prior clustering. `Totem` generates the pseudotime and finds a directed trajectory of the MST using `Slingshot` and can export the trajectory as R objects that are compatible with downstream analysis tools.


## Installation

The easiest way to install Totem is to use devtools R package. We are currently working on to publish the package on Bioconductor.

```R
devtools::install_github("elolab/Totem")
```

## Manual

The HTML version of the R package vignette can be accessed   [here](https://htmlpreview.github.io/?https://github.com/elolab/Totem-benchmarking/blob/main/Totem.html).



## Citation

[1] Johannes Smolander. Sini Junttila. Laura L. Elo. Totem: A user-friendly tool for inferring tree-shaped trajectories from single-cell data. Submitted to *bioRxiv*. 

[2] Street, K., Risso, D., Fletcher, R. et al. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC Genomics 19, 477 (2018). [https://doi.org/10.1186/s12864-018-4772-0](https://doi.org/10.1186/s12864-018-4772-0)

