# Totem

A user-friendly tool for inferring tree-shaped trajectories from single cell data. 

## Description

`Totem` package is a tool designed to facilitate inference of tree-shaped trajectories from single-cell
data [1]. These include linear, bifurcating, multifurcating, and other tree-shaped trajectories. Totem uses Slingshot [2] to find linear curves through the lineages of the tree. 


## Installation

The easiest way to install Totem is to use either `remotes` or `devtools` R package with the `install_github` command. The `dynplot` and  `dynfeature` R packages are not available on CRAN, and they must be installed from GitHub.

```R
remotes::install_github("dynverse/dynfeature")
remotes::install_github("dynverse/dynplot")

remotes::install_github("elolab/Totem")
```

## Manual

The HTML version of the R package vignette can be accessed   [here](https://htmlpreview.github.io/?https://github.com/elolab/Totem-benchmarking/blob/main/Totem.html).


## References

[1] Johannes Smolander, Sini Junttila, Laura L Elo, Cell-connectivity-guided trajectory inference from single-cell data, Bioinformatics, Volume 39, Issue 9, September 2023, btad515, https://doi.org/10.1093/bioinformatics/btad515 

[2] Street, K., Risso, D., Fletcher, R. et al. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC Genomics 19, 477 (2018). [https://doi.org/10.1186/s12864-018-4772-0](https://doi.org/10.1186/s12864-018-4772-0)

