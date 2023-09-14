# Totem

A user-friendly tool for inferring tree-shaped trajectories from single cell data. 

## Description

`Totem` package is a tool designed to facilitate inference of tree-shaped trajectories from single-cell
data [1]. These include linear, bifurcating, multifurcating, and other tree-shaped trajectories. Totem uses Slingshot [2] to find linear curves through the lineages of the tree. 


## Installation

The easiest way to install Totem is to use remotes or devtools R package. The `dynplot` R package is no longer available in CRAN, and it must be installed from GitHub.

```R
remotes::install_github("dynverse/dynplot")
remotes::install_github("elolab/Totem")
```

## Manual

The HTML version of the R package vignette can be accessed   [here](https://htmlpreview.github.io/?https://github.com/elolab/Totem-benchmarking/blob/main/Totem.html).


## References

[1] Johannes Smolander. Sini Junttila. Laura L. Elo. Totem: A user-friendly tool for inferring tree-shaped trajectories from single-cell data. Submitted to *bioRxiv*. 

[2] Street, K., Risso, D., Fletcher, R. et al. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC Genomics 19, 477 (2018). [https://doi.org/10.1186/s12864-018-4772-0](https://doi.org/10.1186/s12864-018-4772-0)

