## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE------------------
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5)

## ----eval = FALSE-------------------------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("Totem")

## -----------------------------------------------------------------------------
set.seed(1234)

suppressMessages(library(Totem))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(cowplot))
suppressMessages(library(dyndimred))
suppressMessages(library(S4Vectors))

data(binary_tree_1)

# The dataset was normalized using the LogNormalize method from the Seurat R package.
sce <- SingleCellExperiment(assays = list(counts = t(binary_tree_1$raw_data),
                                          logcounts = t(binary_tree_1$normalized_data)))
sce <- PrepareTotem(sce)

## -----------------------------------------------------------------------------

print(dim(sce))

sce <- RunDimRed(object = sce,
                 dim.red.method = "lmds",
                 dim.red.features = NULL,
                 dim.reduction.parameter.list = list(ndim=5))


## -----------------------------------------------------------------------------
# extract the generated embedding and set it again into the same slot
own_dim_red <- reducedDim(sce)
reducedDim(sce,type = "lmds") <- own_dim_red # use any name


## -----------------------------------------------------------------------------

dim_red <- dimred_mds(binary_tree_1$normalized_data,ndim=2)


## -----------------------------------------------------------------------------

sce <- RunClustering(sce,
                         k.range = 3:20,
                         min.cluster.size = 5,
                         N.clusterings=10000)


## -----------------------------------------------------------------------------

VizCellConnectivity(sce,custom.dim.red = dim_red)


## -----------------------------------------------------------------------------

sce <- SelectClusterings(sce,selection.method = 1,
                       selection.N.models = 10,
                       selection.stratified=FALSE,
                       prior.clustering = NULL)


## ---- fig.height=12, fig.width=12, out.width = '100%',fig.align = "center"----

ReturnTrajNames(sce)

VizMST(sce,clustering.names = ReturnTrajNames(sce),custom.dim.red = dim_red)


## -----------------------------------------------------------------------------

ReturnMSTNetwork(sce,clustering.name = "13.127")


## -----------------------------------------------------------------------------

sce <- RunSmoothing(sce)


## ---- fig.height=12, fig.width=12, out.width = '100%',fig.align = "center"----

ReturnTrajNames(sce)

VizSmoothedTraj(sce,
                traj.names = ReturnTrajNames(sce),
                custom.dim.red = dim_red,plot.pseudotime = FALSE)


## -----------------------------------------------------------------------------

ReturnSmoothedTrajNetwork(sce,clustering.name = "13.127")


## -----------------------------------------------------------------------------
sce <- ChangeTrajRoot(sce,traj.name="13.127",root.cluster=12)

## ---- fig.height=4, fig.width=4, out.width = '100%',fig.align = "center"------


VizSmoothedTraj(sce,
                traj.names = "13.127",
                custom.dim.red = dim_red,plot.pseudotime = FALSE)

VizSmoothedTraj(sce,
                traj.names = "13.127",
                custom.dim.red = dim_red,plot.pseudotime = TRUE)



## -----------------------------------------------------------------------------

ReturnSmoothedTrajNetwork(sce,clustering.name = "13.127")


## -----------------------------------------------------------------------------
sce <- MergeClusters(sce,clustering.name="13.127",
                     merged.clusters = c(11,3),
                     merged.name = 3)

## -----------------------------------------------------------------------------
VizMST(sce,clustering.names = "13.127",custom.dim.red = dim_red)

## ---- fig.height=4, fig.width=4, out.width = '100%',fig.align = "center"------

VizFeatureExpression(sce,traj.name = "13.127",
                     feature.names = "G112",
                     custom.dim.red = dim_red,
                     plot.trajectory = TRUE)


## ---- fig.height=4, fig.width=4, out.width = '100%',fig.align = "center"------

VizClustering(sce,clustering = binary_tree_1$clustering,custom.dim.red = dim_red)


## -----------------------------------------------------------------------------

dynwrap_object <- ReturnDynwrapObject(sce,traj.name="13.127")
slingshot_object <- ReturnSlingshotObject(sce,traj.name="13.127")



## -----------------------------------------------------------------------------

  milestone_network <- dynwrap_object$milestone_network
  from_ <- milestone_network$from
  to_ <- milestone_network$to
  milestone_network$from <- to_
  milestone_network$to <- from_
  dynwrap_object$milestone_network <- milestone_network
  
  metadata(sce)$totem$dynwrap_trajectory[["13.127"]] <- dynwrap_object



## ---- fig.height=4, fig.width=4, out.width = '100%',fig.align = "center"------


VizSmoothedTraj(sce,
                traj.names = "13.127",
                custom.dim.red = dim_red,plot.pseudotime = FALSE)

VizSmoothedTraj(sce,
                traj.names = "13.127",
                custom.dim.red = dim_red,plot.pseudotime = TRUE)



## ---- eval = TRUE-------------------------------------------------------------
sessionInfo()

