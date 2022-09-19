




#' @title Run Totem
#'
#' @description
#' The function wraps all steps in Totem from dimensionality reduction 
#' to Slingshot smoothing.
#'
#' @param object an object of class \code{SingleCellExperiment}.
#' @param dim.red.method a character denoting which dimensionality 
#' reduction method to use. The options are 
#' `lmds` (Landmark Multi-Dimensional Scaling), 
#' `mds` (Multi-Dimensional Scaling),  #' tsne, umap, pca, 
#' `tsne` (t-distributed stochastic neighbor embedding), 
#' `umap` (Uniform Manifold Approximation and Projection), 
#' `custom` if the user wants to use a custom embedding. The custom
#' embedding can be provided by 
#' dim.reduction.par.list = list(dim_red=X), where `X` is the embedding.
#' @param dim.red.features a character vector specifying the features
#' that are used as input in the dimensionality reduction. Useful in situations
#' where the user wants to use the full gene expression matrix to visualize
#' gene expression but only some of the features in the dimensionality reduction
#' @param dim.reduction.par.list a list specifying the parameters to use
#' in the dimensionality reduction functions of dyndimred R package. 
#' For example, list(ndim=5) can be used to specify the number of 
#' dimensions for lmds and mds.
#' @param k.range a numeric vector specifying 
#' the range of the number of clusters. E.g. 3:20 would generate clusterings
#' for which the number of clusters ranges from 3 to 20.
#' @param min.cluster.size a positive integer specifying 
#' the minimum number of cells the clusters can have. Clusterings that don't
#' meet this requirement are filtered out. Note that this should be generally
#' greater or equal to the number of dimensions 
#' in the dimensionally reduced data. Otherwise the Slingshot distance metric
#' can produce errors.
#' @param N.clusterings a positive integer denoting the initial number of 
#' clusterings that are generated. 
#' After filtering out those with too small clusters the number will go down.
#' @param selection.method a positive integer specifying 
#' the clustering selection method.
#' `1` will use the VRC with the cell connectivity measure as input.
#' `2` will use the VRC with the dimensionally reduced data as input.
#' `3` will use the mean of `1` and `2`.
#' `4` will use the mean of `1` and the ARI between a clustering and 
#' the `prior.clustering`. ARI is the adjusted Rand index that calculates the
#' similarity of two clusterings. `prior.clustering` parameter specifies the
#' prior clustering that the user wants to use to direct the clustering search
#' towards clusterings that are more similar with `prior.clustering`.
#' `5` will use the mean of `2` and the ARI between a clustering and 
#' the `prior.clustering`. ARI is the adjusted Rand index that calculates the
#' similarity of two clusterings. `prior.clustering` parameter specifies the
#' prior clustering that the user wants to use to direct the clustering search
#' towards clusterings that are more similar with `prior.clustering`.
#' @param selection.N.models a positive integer specifying the number of
#' top-ranking clusterings that are selected. If the value is a a numeric 
#' vector, this selects a set of clusterings of which position in the ranked
#' list matches the integers (e.g., 2:4 would select 
#' the clusterings ranked 2nd,3rd and 4th).
#' @param selection.k.range an integer vector that specifies 
#' the numbers of clusters that are allowed in the selection. E.g.
#' c(10,11) would select clusterings that have 10 or 11 clusters.
#' @param selection.stratified a logical variable (TRUE or FALSE) that can be
#' used to use the stratified selection method. If activated, 
#' the number of clusters will vary evenly in the selected clusterings. Can
#' be useful if the top results produce trajectories that are not satisfactory,
#' and the user wants to try very different trajectories.
#' @param prior.clustering a vector that includes a prior clustering. When used
#' in conjunction with `selection.method=4` or `selection.method=5`, it
#' can be used to select clusterings that are more similar with
#' the prior clustering.
#' @param slingshot.par.list a list of the Slingshot parameters.
#' 
#' @name RunTotem
#' @return an object of class \code{SingleCellExperiment}
#'
#' @keywords Totem trajectory inference wrapper
#'
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' 
#' sce <- SingleCellExperiment(assays = list(
#' counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- RunTotem(sce,N.clusterings=200) #use default N.clusterings=10000

RunTotem.SingleCellExperiment <- function(object,
                                          dim.red.method,
                                          dim.red.features,
                                          dim.reduction.par.list,
                                          k.range,
                                          min.cluster.size,
                                          N.clusterings,
                                          selection.method,
                                          selection.N.models,
                                          selection.k.range,
                                          selection.stratified,
                                          prior.clustering,
                                          slingshot.par.list) {
  # Prepare Totem
  object <- PrepareTotem(object)
  
  cat("1. Run dimensionality reduction...")
  
  # Run dimensionality reduction
  object <- RunDimRed(object,dim.red.method,
                      dim.red.features,dim.reduction.par.list)
  
  message("READY!")
  cat("2. Generate a set of clustering results...")
  
  # Run k-means clustering
  object <- RunClustering(object,
                              k.range=k.range,
                              min.cluster.size=min.cluster.size,
                              N.clusterings=N.clusterings)
  
  message("READY!")
  cat("3. Select clusterings...")
  
  object <- SelectClusterings(object,selection.method=selection.method,
                                  selection.N.models=selection.N.models,
                                  selection.k.range=selection.k.range,
                                  selection.stratified=selection.stratified,
                                  prior.clustering=prior.clustering)
  
  message("READY!")
  cat("4. Run the lineage smoothing...")
  
  # Run iterative tree shrinkage
  object <- RunSmoothing(object,slingshot.par.list)
  
  message("READY!")
  return(object)
}

#' @rdname RunTotem
#' @aliases RunTotem
setMethod("RunTotem", signature(object = "SingleCellExperiment"),
          RunTotem.SingleCellExperiment)









#' @title Prepare \code{SingleCellExperiment} object for \code{Totem} analysis
#'
#' @description
#' The function prepares the \code{SingleCellExperiment} object for
#' \code{Totem} analysis. The only required input is an object of class
#' \code{SingleCellExperiment} with data stored at least in the \code{logcounts}
#'  and \code{counts} slots.
#'
#' @param object an object of class \code{SingleCellExperiment}
#'
#' @name PrepareTotem
#'
#' @return an object of class \code{SingleCellExperiment}
#'
#' @keywords prepare totem raw normalized data
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment counts counts<- logcounts logcounts<-
#' @importFrom methods is
#'
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#'
PrepareTotem.SingleCellExperiment <- function(object) {
  
  message("")
  message("The data matrix has ",length(rownames(object)), " features.")
  message("The data matrix has ",length(colnames(object)), " cells.")
  message("Please note that the filtering of poor-quality cells needs to be 
          done before Totem analysis. In addition, the cells need to be part of 
          a single, connected, tree-shaped trajectory 
          (not multiple disconnected trajectories.")
  message("")
  
  # Remove duplicate features
  if (sum(duplicated(rownames(object))) != 0) {
    features_before <- length(rownames(object))
    object <- object[!duplicated(rownames(object)), ]
    features_after <- length(rownames(object))
    message("")
    message("data in SingleCellExperiment object contained duplicate ",
                  " features. ", features_before - features_after,
                  "/", features_before, " were filtered out.")
    message("")
  }
  
  # Filter genes that are not expressed in any of the cells
  genes_before_filtering <- nrow(object)
  non_expressing_genes <- rownames(object)[rowSums(logcounts(object)) != 0]
  object <- object[non_expressing_genes,]
  genes_after_filtering <- nrow(object)
  message("")
  message(genes_after_filtering,"/",genes_before_filtering,
                " features remain after filtering genes with only zero values.")
  message("")
  
  # Create a place into `metadata`` slot for the data from Totem
  metadata(object)$totem <- list()
  
  return(object)
}

#' @rdname PrepareTotem
#' @aliases PrepareTotem
setMethod("PrepareTotem", signature(object = "SingleCellExperiment"),
          PrepareTotem.SingleCellExperiment)



#' @title Run dimensionality reduction
#'
#' @description
#' The function performs dimensionality reduction using 
#' \code{dyndimred} package methods.
#' @param object an object of class \code{SingleCellExperiment}.
#' @param dim.red.method a character denoting which dimensionality 
#' reduction method to use. The options are 
#' `lmds` (Landmark Multi-Dimensional Scaling), 
#' `mds` (Multi-Dimensional Scaling),  #' tsne, umap, pca
#' `tsne` (t-distributed stochastic neighbor embedding), 
#' `umap` (Uniform Manifold Approximation and Projection), 
#' `custom` if the user wants to use a custom embedding. The custom
#' embedding can be provided by 
#' dim.reduction.par.list = list(dim_red=X), where `X` is the embedding.
#' @param dim.red.features a character vector specifying the features
#' that are used as input in the dimensionality reduction. Useful in situations
#' where the user wants to use the full gene expression matrix to visualize
#' gene expression but only some of the features in the dimensionality reduction
#' @param dim.reduction.par.list a list specifying the parameters to use
#' in the dimensionality reduction functions of dyndimred R package. 
#' For example, list(ndim=5) can be used to specify the number of 
#' dimensions for lmds and mds.
#' @name RunDimRed
#'
#' @return an object of class \code{SingleCellExperiment}
#'
#' @keywords dimensionality reduction features
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
#' @importFrom dyndimred dimred_pca dimred_landmark_mds dimred_mds
#' @importFrom dyndimred dimred_umap dimred_tsne
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))

RunDimRed.SingleCellExperiment <- function(object,dim.red.method,
                                           dim.red.features,
                                           dim.reduction.par.list) {
  
  X <- Matrix::t(logcounts(object))
  
  
  if (!is.null(dim.red.features)) {
    message("using features in the 'dim.red.features' argument")
    
    X <- X[dim.red.features,]
  }
  
  message("using '",dim.red.method,"' as the dimensionality reduction method")
  
  if (dim.red.method=="tsne") {
    X <- dimred_tsne(X,ndim = dim.reduction.par.list[["ndim"]],
                     initial_dims=dim.reduction.par.list[["initial_dims"]])
    reducedDim(object,type = "tsne") <- X
  }
  else if (dim.red.method=="umap") {
    X <- dimred_umap(X,ndim = dim.reduction.par.list[["ndim"]],
                     pca_components=dim.reduction.par.list[["pca_components"]])
    reducedDim(object,type = "umap") <- X
  }
  else if (dim.red.method=="pca") {
    X <- dimred_pca(X,ndim = dim.reduction.par.list[["ndim"]])
    reducedDim(object,type = "pca") <- X
  }
  else if (dim.red.method=="lmds") {
    X <- dimred_landmark_mds(X,ndim = dim.reduction.par.list[["ndim"]])
    reducedDim(object,type = "lmds") <- X
    
  } 
  else if (dim.red.method=="mds") {
    X <- dimred_mds(X,ndim = dim.reduction.par.list[["ndim"]])
    reducedDim(object,type = "mds") <- X
  } 
  else if(dim.red.method=="custom") {
    message("Using dim.reduction.par.list[['dim_red']] ",
            "as the dimensionally reduced data")
    X <- dim.reduction.par.list[["dim_red"]]
    reducedDim(object,type = "custom") <- X
    
  } 
  else {
    message("No dimensionality reduction performed!")
  }
  return(object)
}

#' @rdname RunDimRed
#' @aliases RunDimRed
setMethod("RunDimRed", signature(object = "SingleCellExperiment"),
          RunDimRed.SingleCellExperiment)





#' @title Run k-medoids clustering
#'
#' @description
#' The function runs the k-medoids clustering, filters out clusterings that
#' have too small clusters, and calculates the cell connectivity measure.
#' @param object an object of class \code{SingleCellExperiment}
#' @param k.range a numeric vector specifying 
#' the range of the number of clusters. E.g. 3:20 would generate clusterings
#' for which the number of clusters ranges from 3 to 20.
#' @param min.cluster.size a positive integer specifying 
#' the minimum number of cells the clusters can have. Clusterings that don't
#' meet this requirement are filtered out. Note that this should be generally
#' greater or equal to the number of dimensions 
#' in the dimensionally reduced data. Otherwise the Slingshot distance metric
#' can produce errors.
#' @param N.clusterings a positive integer denoting the initial number of 
#' clustering reuslts. After filtering out those with too small clusters the 
#' number will go down.
#' @param verbose a logical for denoting whether the index of the clustering
#' that is generated should be printed for progress following purposes.
#' @name RunClustering
#'
#' @return an object of class \code{SingleCellExperiment}
#'
#' @keywords clustering clara k-medoids
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#' @importFrom fpc calinhara
#'
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000

RunClustering.SingleCellExperiment <- function(object,
                                               k.range,
                                               min.cluster.size,
                                               N.clusterings,
                                               verbose) {
  
  if (identical(reducedDimNames(object), character(0))) {
    message("No dimensionally reduced data available.",
            " Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- as.matrix(Matrix::t(X))
  } else {
    X <- reducedDim(object)
    message("using ", reducedDimNames(object), " as the dim.reduced data")
  }
  ks <- sample(k.range,size = N.clusterings,replace = TRUE)
  ks <- make.unique(as.character(ks))
  
  clusterings <- RunCLARA(X,ks,min.cluster.size,verbose)
  ks <- names(clusterings)
  msts <- GenerateMST(X,clusterings)
  connectivities <- Connectivity(clusterings,msts)
  cellcon <- apply(do.call(cbind,connectivities),1,mean)
  
  vrcs_conn <- unlist(lapply(clusterings,function(x) calinhara(cellcon,x)))
  vrcs_conn <- vrcs_conn/max(vrcs_conn)
  
  vrcs_X <- lapply(ks,function(x) calinhara(x=X,clustering=clusterings[[x]]))
  vrcs_X <- unlist(vrcs_X)
  names(vrcs_X) <- ks
  vrcs_X <- vrcs_X/max(vrcs_X)
  
  mean_vrc_cellconn_vrcs_X <- (vrcs_conn+vrcs_X)/2
  
  df <- cbind(vrcs_conn,vrcs_X,mean_vrc_cellconn_vrcs_X,ks)
  rownames(df) <- ks
  colnames(df) <- c("vrc_cellconn","vrc_dimred",
                    "mean_vrc_cellconn_vrc_dimred","clustering_name")
  
  df <- as.data.frame(df)
  
  df$vrc_cellconn <- as.numeric(df$vrc_cellconn)
  df$vrc_dimred <- as.numeric(df$vrc_dimred)
  df$mean_vrc_cellconn_vrc_dimred <- as.numeric(df$mean_vrc_cellconn_vrc_dimred)

  metadata(object)$totem$cell.connectivity <- cellcon
  metadata(object)$totem$all.clustering <- clusterings
  metadata(object)$totem$all.clustering.scores <- df
  
  return(object)
}

#' @rdname RunClustering
#' @aliases RunClustering
setMethod("RunClustering", signature(object = "SingleCellExperiment"),
          RunClustering.SingleCellExperiment)





#' @title Perform clustering selection
#'
#' @description
#' The function selects clusterings that are used in the 
#' trajectory smoothing with Slingshot. The user can choose from
#' different criteria to perform the selection.
#' @param object an object of class \code{SingleCellExperiment}
#' @param selection.method a positive integer specifying 
#' the clustering selection method.
#' \code{1} will use the Variance Ratio Criterion (Calinski-Harabasz score) 
#' VRC with the cell connectivity measure as input.
#' \code{2} will use the VRC with the dimensionally reduced data as input.
#' \code{3} will use the mean of \code{1} and \code{2}.
#' \code{4} will use the mean of \code{1} and the ARI between a clustering and 
#' the \code{prior.clustering}. ARI is the adjusted Rand index that calculates 
#' the similarity of two clusterings. \code{prior.clustering} parameter 
#' specifies the prior clustering that the user wants to use to direct the 
#' clustering search towards clusterings that are more similar with 
#' \code{prior.clustering}. \code{5} will use the mean of \code{2} and 
#' the ARI between a clustering and 
#' the \code{prior.clustering}. ARI is the adjusted Rand index that measures the
#' similarity of two clusterings.\code{prior.clustering} parameter specifies the
#' prior clustering that the user wants to use to direct the clustering search
#' towards clusterings that are more similar with \code{prior.clustering}.
#' @param selection.N.models a positive integer specifying the number of
#' top-ranking clusterings that are selected. If the value is a a numeric 
#' vector, this selects a set of clusterings of which position in the ranked
#' list matches the integers (e.g., \code{2:4} would select 
#' the clusterings ranked 2nd,3rd and 4th).
#' @param selection.k.range an integer vector that specifies 
#' the range of the number of clusters that are allowed in the selection. E.g.
#' \code{c(10,11)} would select clusterings that have 10 or 11 clusters.
#' @param selection.stratified a logical variable (TRUE or FALSE) that can be
#' used to activate the stratified selection method. If activated, 
#' the number of clusters will vary evenly in the selected clusterings. Can
#' be useful if the top results include trajectories that are not satisfactory,
#' and the user wants to try diverse trajectories.
#' @param prior.clustering a vector that includes a prior clustering. When used
#' in conjunction with \code{selection.method=4} or \code{selection.method=5},
#' it can be used to select clusterings that are more similar with
#' the prior clustering.
#' 
#' @name SelectClusterings
#'
#' @return an object of class \code{SingleCellExperiment}
#'
#' @keywords clustering selection ARI
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#' @importFrom aricode ARI
#'
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' sce <- SelectClusterings(sce)

SelectClusterings.SingleCellExperiment <- function(object,
                                                       selection.method,
                                                       selection.N.models,
                                                       selection.k.range,
                                                       selection.stratified,
                                                       prior.clustering) {
  clusterings <- metadata(object)$totem$all.clustering
  all.clust.scores <- metadata(object)$totem$all.clustering.scores
  ks <- names(clusterings)
  if (!is.null(selection.k.range)) {
    ks <- ks[as.numeric(gsub("\\..*","",ks)) %in% selection.k.range]
    clusterings <- clusterings[ks]
    all.clust.scores <- all.clust.scores[ks,]
  }
  if (is.null(all.clust.scores$ARI) && selection.method > 3) {
    ari <- unlist(lapply(metadata(object)$totem$all.clustering,
                         function(x) ARI(x,prior.clustering)))
    metadata(object)$totem$all.clustering.scores$ARI <- ari
    all.clust.scores$ARI <- ari[ks]
  }
  if (selection.method==1) {
    ave_score <- all.clust.scores$vrc_cellconn
  } else if (selection.method==2) {
    ave_score <- all.clust.scores$vrc_dimred
  } else if (selection.method==3) {
    ave_score <- all.clust.scores$mean_vrc_cellconn_vrc_dimred
  } else if (selection.method==4) {
    ave_score <- (all.clust.scores$vrc_cellconn+all.clust.scores$ARI)/2
  } else if (selection.method==5) {
    ave_score <- (all.clust.scores$vrc_dimred+all.clust.scores$ARI)/2
  }
  names(ave_score) <- ks
  k.range <- names(table(as.numeric(gsub("\\..*","",ks))))
  ave_score <- sort(ave_score,decreasing = TRUE)
  ave_score_ks <- names(ave_score)
  if (!selection.stratified) {
    if (length(selection.N.models)==1) {
      clusterings_ <- clusterings[ave_score_ks][seq_len(selection.N.models)]
    } else {
      clusterings_ <- clusterings[ave_score_ks][selection.N.models]
    }
  } else {
    ave_score_ks_split <- split(ave_score_ks,
                                f = as.numeric(gsub("\\..*","",ave_score_ks)))
    k_included <- c()
    k <- min(as.numeric(names(ave_score_ks_split)))
    while (length(k_included) < max(selection.N.models)) {
      k_selected <- unlist(lapply(ave_score_ks_split,
                                  function(x) {setdiff(x,k_included)[1]}))
      k_selected <- k_selected[order(ave_score[k_selected],decreasing = TRUE)]
      k_included <- c(k_included,k_selected)
    } 
    if (length(selection.N.models) == 1) {
      k_included <- k_included[seq(1,selection.N.models)]
    } else {
      k_included <- k_included[selection.N.models]
    }
    clusterings_ <- clusterings[k_included]
  }
  metadata(object)$totem$selected.clustering <- clusterings_
  return(object)
}

#' @rdname SelectClusterings
#' @aliases SelectClusterings
setMethod("SelectClusterings", signature(object = "SingleCellExperiment"),
          SelectClusterings.SingleCellExperiment)


#' @title Run Slingshot smoothing
#'
#' @description
#' The function runs the Slingshot smoothing. The Slingshot parameters can also
#' be adjusted. The Slignshot smoothing is performed for the selected 
#' clusterings that were selected using \code{SelectClusterings} function.
#'
#' @param object an object of class \code{SingleCellExperiment}.
#' @param slingshot.par.list a list of the Slingshot parameters.
#'
#' @name RunSmoothing
#' @return an object of class \code{SingleCellExperiment}
#'
#' @keywords slingshot smoothing
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' 
#'
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) #Use default N.clusterings=10000
#' sce <- SelectClusterings(sce)
#' sce <- RunSmoothing(sce)

RunSmoothing.SingleCellExperiment <- function(object,slingshot.par.list) {
  
  clusterings <- metadata(object)$totem$selected.clustering
  
  if (identical(reducedDimNames(object), character(0)))
  {
    message("No dimensionality reduced data available. ",
            "Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- as.matrix(Matrix::t(X))
  } else
  {
    X <- reducedDim(object)
  }
  
  dynwrap_objects <- list()
  slingshot_objects <- list()
  for (kkk in names(clusterings))
  {
    clustering <- clusterings[[kkk]]
    
    traj <- try(RunSlingshotSmoothing(rd=X,
                                      labels=clustering,
                                      start.clus=NULL,
                                      slingshot.par.list=slingshot.par.list))
    if (is(traj,"try-error"))
    {
      traj <- NULL
    }
    
    dynwrap_objects[[kkk]] <- traj$dynwrap
    slingshot_objects[[kkk]] <- traj$slingshot
    
  }
  metadata(object)$totem$dynwrap_trajectory <- dynwrap_objects
  metadata(object)$totem$slingshot_trajectory <- slingshot_objects
  
  
  return(object)
}

#' @rdname RunSmoothing
#' @aliases RunSmoothing
setMethod("RunSmoothing", signature(object = "SingleCellExperiment"),
          RunSmoothing.SingleCellExperiment)






#' @title Visualize a Minimum Spanning Tree
#'
#' @description
#' The function visualizes a Minimum Spanning Tree on an embedding. The
#' MST is generated based on the clusterings of which names are
#' given through the \code{clustering.names} parameter. The distance method
#' is the Slingshot distance method, which can not be currently changed.
#' @param object an object of class \code{SingleCellExperiment}.
#' @param clustering.names a character vector that specifies the names of the
#' MSTs (clusterings) to be visualized.
#' @param dim.red.type a character that denotes the name of the 
#' dimensionality reduction method stored in the \code{SingleCellExperiment}
#' object.
#' @param viz.dim.red a matrix that specifies a custom two-dimensional
#' embedding for visualization.
#'
#' @name VizMST
#'
#' @return a \code{ggplot2} object
#'
#' @keywords MST visualization
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#' @importFrom cowplot plot_grid
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' sce <- SelectClusterings(sce)
#' sce <- RunSmoothing(sce)
#' dim_red <- dyndimred::dimred_mds(t(logcounts(sce)),ndim=2)
#' clustering_name <- ReturnTrajNames(sce)[1]
#' VizMST(sce,viz.dim.red=dim_red,clustering.names=clustering_name)

VizMST.SingleCellExperiment <- function(object,
                                        clustering.names,
                                        dim.red.type,
                                        viz.dim.red) {
  
  if (identical(reducedDimNames(object), character(0))) {
    message("No dimensionality reduced data available. ",
            "Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- as.matrix(Matrix::t(X))
  } else {
    X <- reducedDim(object)
  }
  if (!is.null(dim.red.type)) {
    X <- reducedDim(object,type = dim.red.type)
  }
  
  plot_list <- list()
  for (clustering_name in clustering.names)
  {
    
    clustering <- metadata(object)$totem$all.clustering[[clustering_name]]
    
    p <- VizOneMST(viz.dim.red=viz.dim.red,
                   X=X,
                   clustering=clustering,
                   clustering.name=clustering_name)
    
    plot_list[[clustering_name]] <- p
  }
  
  p <- plot_grid(plotlist = plot_list)
  return(p)
  
}

#' @rdname VizMST
#' @aliases VizMST
setMethod("VizMST", signature(object = "SingleCellExperiment"),
          VizMST.SingleCellExperiment)





#' @title Visualize a user-provided clustering
#'
#' @description
#' The function can be used to visualize a user-provided clustering.
#' The MST of the clustering can also be visualized, which is generated
#' using the Slingshot distance method.
#' @param object an object of class \code{SingleCellExperiment}.
#' @param clustering a vector that specifies the clustering.
#' @param dim.red.type a character that denotes the name of the 
#' dimensionality reduction method stored in the \code{SingleCellExperiment}
#' object.
#' @param viz.dim.red a matrix that specifies a custom two-dimensional
#' embedding for visualization.
#' @param plot.mst a logical denoting whether the MST should
#' be visualized. 
#' 
#' @name VizClustering
#'
#' @return a \code{ggplot2} object
#'
#' @keywords visualization embedding two-dimensional clustering
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#' @importFrom dynplot plot_dimred
#' @importFrom ggplot2 ggtitle theme element_text
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' sce <- SelectClusterings(sce)
#' sce <- RunSmoothing(sce)
#' dim_red <- dyndimred::dimred_mds(t(logcounts(sce)),ndim=2)
#' clustering_name <- ReturnTrajNames(sce)[1]
#' clustering <- ReturnClustering(sce,clustering_name)
#' VizClustering(sce,viz.dim.red=dim_red,clustering)

VizClustering.SingleCellExperiment <- function(object,
                                               clustering,
                                               dim.red.type,
                                               viz.dim.red,
                                               plot.mst) {
  
  if (identical(reducedDimNames(object), character(0)))
  {
    message("No dimensionality reduced data available. ",
            "Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- as.matrix(Matrix::t(X))
  } else
  {
    X <- reducedDim(object)
  }
  
  if (!is.null(dim.red.type))
  {
    X <- reducedDim(object,type = dim.red.type)
  }
  
  p <- VizOneMST(viz.dim.red=viz.dim.red,
                 X=X,
                 clustering=clustering,
                 clustering.name="")
  
  return(p)
  
}

#' @rdname VizClustering
#' @aliases VizClustering
setMethod("VizClustering", signature(object = "SingleCellExperiment"),
          VizClustering.SingleCellExperiment)





#' @title Visualize a smoothed trajectory
#'
#' @description
#' The function visualizes a Slingshot-smoothed trajectory.
#' @param object an object of class \code{SingleCellExperiment}.
#' @param traj.names a character vector of the trajectory (clustering) names.
#' @param dim.red.type a character that specifies the name of the 
#' dimensionality reduction method saved in the \code{SingleCellExperiment} 
#' object.
#' @param viz.dim.red a matrix for the user-provided, custom dimensionality
#' reduction.
#' @param plot.pseudotime a logical for whether the pseudotime of 
#' the trajectories should be plotted.
#' 
#' @name VizSmoothedTraj
#'
#' @return a \code{ggplot2} object
#'
#' @keywords visualization dimensionality reduction embedding smooth trajectory
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#' @importFrom dynplot plot_dimred
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggtitle theme element_text
#'
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' sce <- SelectClusterings(sce)
#' sce <- RunSmoothing(sce)
#' dim_red <- dyndimred::dimred_mds(t(logcounts(sce)),ndim=2)
#' clustering_name <- ReturnTrajNames(sce)[1]
#' VizSmoothedTraj(sce,viz.dim.red=dim_red,traj.names=clustering_name)

VizSmoothedTraj.SingleCellExperiment <- function(object,
                                                 traj.names,
                                                 dim.red.type,
                                                 viz.dim.red,
                                                 plot.pseudotime) {
  
  if (identical(reducedDimNames(object), character(0)))
  {
    message("No dimensionality reduced data available. ",
            "Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- as.matrix(Matrix::t(X))
  } else
  {
    X <- reducedDim(object)
  }
  if (!is.null(dim.red.type))
  {
    X <- reducedDim(object,type = dim.red.type)
  }
  plot_list <- list()
  for (traj_name in traj.names) {
    model <- metadata(object)$totem$dynwrap_trajectory[[traj_name]]
    if (plot.pseudotime) {
      plot_list[[traj_name]] <- plot_dimred(model,
                                            color_cells = "pseudotime",
                                            dimred = viz.dim.red) + 
        ggtitle(traj_name)
    } else {
      plot_list[[traj_name]] <- plot_dimred(model,
                                            label_milestones = TRUE,
                                            dimred = viz.dim.red) + 
        ggtitle(traj_name)
    }
  }
  p <- cowplot::plot_grid(plotlist = plot_list)
  return(p)
}

#' @rdname VizSmoothedTraj
#' @aliases VizSmoothedTraj
setMethod("VizSmoothedTraj", signature(object = "SingleCellExperiment"),
          VizSmoothedTraj.SingleCellExperiment)





#' @title Visualize feature expression
#'
#' @description
#' The function visualizes expression of one or multiple features
#' on an embedding .
#' @param object an object of class \code{SingleCellExperiment}.
#' @param traj.name a character for the trajectory (clustering) name.
#' @param feature.names a character vector of the feature names.
#' @param dim.red.type a character that specifies the name of the dimensionality
#'  reduction method saved in the \code{SingleCellExperiment} object.
#' @param viz.dim.red a matrix for the user-provided, custom dimensionality
#' reduction.
#' @param plot.traj a logical for whether the trajectory should 
#' be plotted.
#' 
#' @name VizFeatureExpression
#'
#' @return a \code{ggplot2} object
#'
#' @keywords visualization feature expression dimensionality reduction embedding
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#' @importFrom dynplot plot_dimred
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggtitle theme element_text
#'
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' sce <- SelectClusterings(sce)
#' sce <- RunSmoothing(sce)
#' dim_red <- dyndimred::dimred_mds(t(logcounts(sce)),ndim=2)
#' clustering_name <- ReturnTrajNames(sce)[1]
#' feature_names <- rownames(logcounts(sce))[c(1,2)]
#' VizFeatureExpression(sce,viz.dim.red=dim_red,feature.names=feature_names)

VizFeatureExpression.SingleCellExperiment <- function(object,
                                                      traj.name,
                                                      feature.names,
                                                      dim.red.type,
                                                      viz.dim.red,
                                                      plot.traj) {
  
  if (identical(reducedDimNames(object), character(0))) {
    message("No dimensionality reduced data available. ",
            "Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- as.matrix(Matrix::t(X))
  } else {
    X <- reducedDim(object)
  }
  if (!is.null(dim.red.type)) {
    X <- reducedDim(object,type = dim.red.type)
  }
  exp_matrix <- Matrix::t(logcounts(object))
  plot_list <- list()
  if (is.null(traj.name)) {
    # create random model
    clustering <- metadata(object)$totem$all.clustering[[1]]
    p_traj <- VizOneMST(viz.dim.red=viz.dim.red,
                   X=X,
                   clustering=clustering,
                   clustering.name="",
                   return.traj=TRUE)
    model <- p_traj$traj
    for (feature_name in feature.names) {
      plot_list[[feature_name]] <- plot_dimred(model,
                                               dimred = viz.dim.red,
                                               feature_oi = feature_name,
                                               expression_source = exp_matrix,
                                               plot_trajectory = plot.traj)
    }
  } else {
    model <- metadata(object)$totem$dynwrap_trajectory[[traj.name]]
    for (feature_name in feature.names) {
      plot_list[[feature_name]] <- plot_dimred(model,
                                               dimred = viz.dim.red,
                                               feature_oi = feature_name,
                                               expression_source = exp_matrix,
                                               plot_trajectory = plot.traj) + 
        ggtitle(traj.name)
    }
  }
  return(plot_grid(plotlist = plot_list))
}

#' @rdname VizFeatureExpression
#' @aliases VizFeatureExpression
setMethod("VizFeatureExpression", signature(object = "SingleCellExperiment"),
          VizFeatureExpression.SingleCellExperiment)





#' @title Visualize the cell connectivity measure
#'
#' @description
#' The function visualizes the cell connectivity measure, which is used as 
#' a basis for determining the shape of the trajectory. A low value indicates 
#' that the cell population has, on average, a small number of connections 
#' with other cell types. These are generally the leaf/end nodes of 
#' the trajectory tree. Cell subsets that have a high value are 
#' at the branching points of the trajectory, and they have, on average, more
#' connections compared to the other cell types. The measure is calculated by
#' counting in each MST how many edges the clusters have, divided by the number
#' of clusters in each clustering, and by averaging the ratios of all the MSTs.
#' @param object an object of class \code{SingleCellExperiment}.
#' @param dim.red.type a character that specifies the name of the dimensionality
#'  reduction method saved in the \code{SingleCellExperiment} object.
#' @param viz.dim.red a matrix for the user-provided, custom dimensionality
#' reduction.
#' 
#' @name VizCellConnectivity
#'
#' @return a \code{ggplot2} object
#'
#' @keywords cell connectivity visualization
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#' @importFrom dynplot plot_dimred
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggtitle theme element_text
#' @importFrom ape mst
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) #Use default N.clusterings=10000
#' dim_red <- dyndimred::dimred_mds(t(logcounts(sce)),ndim=2)
#' VizCellConnectivity(sce,viz.dim.red=dim_red)

VizCellConnectivity.SingleCellExperiment <- function(object,
                                                         dim.red.type,
                                                         viz.dim.red) {
  
  if (identical(reducedDimNames(object), character(0))) {
    message("No dimensionality reduced data available. ",
            "Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- as.matrix(Matrix::t(X))
  } else {
    X <- reducedDim(object)
  }
  if (!is.null(dim.red.type)) {
    X <- reducedDim(object,type = dim.red.type)
  }
  cellcon <- metadata(object)$totem$cell.connectivity
  
  expression_matrix <- as.matrix(Matrix::t(logcounts(object)))
  feature_names <- colnames(expression_matrix)
  expression_matrix <- cbind(expression_matrix,cellcon)
  colnames(expression_matrix) <- c(feature_names,"cell_connectivity")
  
  clustering <- metadata(object)$totem$all.clustering[[1]]
  
  p_traj <- VizOneMST(viz.dim.red=viz.dim.red,
                      X=X,
                      clustering=clustering,
                      clustering.name="",
                      return.traj=TRUE)
  model <- p_traj$traj
  
  p <- plot_dimred(model,
                   dimred = viz.dim.red,
                   feature_oi = "cell_connectivity",
                   expression_source = expression_matrix,
                   plot_trajectory = FALSE)
  
  return(p)
  
}

#' @rdname VizCellConnectivity
#' @aliases VizCellConnectivity
setMethod("VizCellConnectivity", signature(object = "SingleCellExperiment"),
          VizCellConnectivity.SingleCellExperiment)





#' @title Return trajectory names
#'
#' @description
#' The function returns the names of the trajectories, 
#' i.e. those that were selected using the \code{SelectClusterings} function.
#' @param object an object of class \code{SingleCellExperiment}
#' @name ReturnTrajNames
#'
#' @return a character vector
#'
#' @keywords clustering names
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' sce <- SelectClusterings(sce,selection.method = 1,selection.N.models = 2)
#' ReturnTrajNames(sce)

ReturnTrajNames.SingleCellExperiment <- function(object) {
  
  clustering_names <- names(metadata(object)$totem$selected.clustering)
  return(clustering_names)
  
}

#' @rdname ReturnTrajNames
#' @aliases ReturnTrajNames
setMethod("ReturnTrajNames", signature(object = "SingleCellExperiment"),
          ReturnTrajNames.SingleCellExperiment)











#' @title Return a clustering
#'
#' @description
#' The function returns a clustering corresponding 
#' to the provided clustering name.
#' @param object an object of class \code{SingleCellExperiment}
#' @param clustering.name a character for the clustering name.
#' 
#' @name ReturnClustering
#'
#' @return a named character vector
#'
#' @keywords clustering
#'
#' @importFrom S4Vectors metadata metadata<-
#'
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' sce <- SelectClusterings(sce,selection.method = 1,selection.N.models = 2)
#' clustering_name <-  ReturnTrajNames(sce)[1]
#' clustering <- ReturnClustering(sce,clustering_name)

ReturnClustering.SingleCellExperiment <- function(object,clustering.name) {
  
  return(metadata(object)$totem$all.clustering[[clustering.name]])
  
}

#' @rdname ReturnClustering
#' @aliases ReturnClustering
setMethod("ReturnClustering", signature(object = "SingleCellExperiment"),
          ReturnClustering.SingleCellExperiment)

















#' @title Return an MST network of a clustering
#'
#' @description
#' The function returns the MST network of a clustering.
#' @param object an object of class \code{SingleCellExperiment}
#' @param clustering.name a character for the clustering name.
#' @name ReturnMSTNetwork
#'
#' @return an object of class "mst"
#'
#' @keywords MST network clustering
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#' @importFrom ape mst
#'
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' clustering_name <-  "3.1"
#' mst_clustering <- ReturnMSTNetwork(sce,clustering_name)

ReturnMSTNetwork.SingleCellExperiment <- function(object,clustering.name) {
  
  if (identical(reducedDimNames(object), character(0)))
  {
    message("No dimensionality reduced data available. ",
            "Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- as.matrix(Matrix::t(X))
  } else
  {
    X <- reducedDim(object)
    message("using ", reducedDimNames(object), 
            " as the dimensionally reduced data")
  }
  
  clustering <- metadata(object)$totem$all.clustering[[clustering.name]]
  
  centroids <- ClusterCentroids(clustering,X)
  mst_out <- mst(SlingshotDistanceMatrix(x=X,centers=centroids,
                                         clusters=clustering))
  
  return(mst_out)
  
}

#' @rdname ReturnMSTNetwork
#' @aliases ReturnMSTNetwork
setMethod("ReturnMSTNetwork", signature(object = "SingleCellExperiment"),
          ReturnMSTNetwork.SingleCellExperiment)




#' @title Return the milestone network of a trajectory
#'
#' @description
#' The function returns the milestone network of a trajectory 
#' smoothed with Slingshot.
#' @param object an object of class \code{SingleCellExperiment}
#' @param clustering.name a character denoting 
#' the name of the trajectory (clustering).
#' @name ReturnSmoothedTrajNetwork
#'
#' @return a data frame
#'
#' @keywords trajectory network
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#'
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) #Use default N.clusterings=10000
#' sce <- SelectClusterings(sce,selection.method = 1,selection.N.models = 2)
#' sce <- RunSmoothing(sce)
#' clustering_name <-  ReturnTrajNames(sce)[1]
#' mst_clustering <- ReturnMSTNetwork(sce,clustering_name)

ReturnSmoothedTrajNetwork.SingleCellExperiment <- function(object,
                                                           clustering.name) {
  
  dynwrap_objects <- metadata(object)$totem$dynwrap_trajectory
  return(as.data.frame(dynwrap_objects[[clustering.name]]$milestone_network))
  
}

#' @rdname ReturnSmoothedTrajNetwork
#' @aliases ReturnSmoothedTrajNetwork
setMethod("ReturnSmoothedTrajNetwork", 
          signature(object = "SingleCellExperiment"),
          ReturnSmoothedTrajNetwork.SingleCellExperiment)












#' @title Change the root cluster of a trajectory
#'
#' @description
#' The function changes the root cluster of a trajectory, which is initially
#' chosen randomly.
#' @param object an object of class \code{SingleCellExperiment}
#' @param traj.name a character for the name of the trajectory
#' @param root.cluster an integer for the name of the new root cluster
#' @name ChangeTrajRoot
#'
#' @return an object of class \code{SingleCellExperiment}
#'
#' @keywords root trajectory start node cluster
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#'
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' sce <- SelectClusterings(sce,selection.method = 1,selection.N.models = 2)
#' sce <- RunSmoothing(sce)
#' clustering_name <-  ReturnTrajNames(sce)[1]
#' sce <- ChangeTrajRoot(sce,traj.name=clustering_name,root.cluster=2)

ChangeTrajRoot.SingleCellExperiment <- function(object,traj.name,root.cluster) {
  
  
  clusterings <- metadata(object)$totem$selected.clustering
  
  if (identical(reducedDimNames(object), character(0))) {
    message("No dimensionality reduced data available. ",
            "Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- as.matrix(Matrix::t(X))
  } else {
    X <- reducedDim(object)
  }
  clustering <- clusterings[[traj.name]]
  
  old_trajectory <- metadata(object)$totem$slingshot_trajectory[[traj.name]]
  slingshot_parameters <- old_trajectory@metadata$slingParams
  
  traj <- try(RunSlingshotSmoothing(X,clustering,root.cluster,
                                    slingshot_parameters))
  if (is(traj,"try-error")) {
    traj <- NULL
  }
  
  metadata(object)$totem$dynwrap_trajectory[[traj.name]] <- traj$dynwrap
  metadata(object)$totem$slingshot_trajectory[[traj.name]] <- traj$slingshot
  
  return(object)
}

#' @rdname ChangeTrajRoot
#' @aliases ChangeTrajRoot
setMethod("ChangeTrajRoot", signature(object = "SingleCellExperiment"),
          ChangeTrajRoot.SingleCellExperiment)









#' @title Merge clusters of a clustering
#'
#' @description
#' The function merges clusters of a clustering. The original clustering will be
#' replaced by the new one, and the only way to get it back is to re-run the 
#' \code{RunClustering} function.
#' @param object an object of class \code{SingleCellExperiment}
#' @param clustering.name a character for the name of the clustering.
#' @param merged.clusters a vector denoting the names of the clusters to merge.
#' @param merged.name a character denoting the name of the merged clusters.
#' By default, the new name will be the name of the cluster that appears first
#' in the clustering.
#' @name MergeClusters
#' @return an object of class \code{SingleCellExperiment}
#'
#' @keywords merge clustering clusters
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom plyr mapvalues
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) #Use default N.clusterings=10000
#' sce <- SelectClusterings(sce,selection.method = 1,selection.N.models = 2)
#' clustering_name <-  ReturnTrajNames(sce)[1]
#' sce <- MergeClusters(sce,clustering.name=clustering_name,
#'                      merged.clusters=c(1,2),merged.name=c(1))

MergeClusters.SingleCellExperiment <- function(object,clustering.name,
                                               merged.clusters,merged.name) {
  
  
  clusterings <- metadata(object)$totem$all.clustering
  
  clustering <- clusterings[[clustering.name]]
  
  if (is.null(merged.name))
  {
    merged.name <- clustering[1]
  }
  
  clustering <- mapvalues(clustering,from=merged.clusters,
                          to = rep(merged.name,length(merged.clusters)))
  metadata(object)$totem$all.clustering[[clustering.name]] <- clustering
  if (clustering.name %in% names(metadata(object)$totem$selected.clustering))
  {
    metadata(object)$totem$selected.clustering[[clustering.name]] <- clustering
  }

  return(object)
  
}

#' @rdname MergeClusters
#' @aliases MergeClusters
setMethod("MergeClusters", signature(object = "SingleCellExperiment"),
          MergeClusters.SingleCellExperiment)













#' @title Return a trajectory as a dynwrap object
#'
#' @description
#' The function returns a Slingshot-smoothed trajectory 
#' as a \code{dynwrap} object
#' @param object an object of class \code{SingleCellExperiment}
#' @param traj.name a character for the trajectory name
#' 
#' @name ReturnDynwrapObject
#'
#' @return an object of class \code{dynwrap}
#'
#' @keywords dynwrap trajectory smooth
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#'
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' sce <- SelectClusterings(sce,selection.method = 1,selection.N.models = 2)
#' sce <- RunSmoothing(sce)
#' clustering_name <-  ReturnTrajNames(sce)[1]
#' dynwrap_object <- ReturnDynwrapObject(sce,clustering_name)

ReturnDynwrapObject.SingleCellExperiment <- function(object,traj.name) {
  
  
  return(metadata(object)$totem$dynwrap_trajectory[[traj.name]])
  
  
  
}

#' @rdname ReturnDynwrapObject
#' @aliases ReturnDynwrapObject
setMethod("ReturnDynwrapObject", signature(object = "SingleCellExperiment"),
          ReturnDynwrapObject.SingleCellExperiment)












#' @title Return a trajectory as a PseudotimeOrdering object
#'
#' @description
#' The function returns a smoothed trajectory as a \code{PseudotimeOrdering} 
#' object, which is the object class of Slingshot.
#' @param object an object of class \code{SingleCellExperiment}
#' @param traj.name a character for the trajectory name
#' 
#' @name ReturnSlingshotObject
#'
#' @return an object of class \code{PseudotimeOrdering}
#'
#' @keywords PseudotimeOrdering Slingshot object trajectory smooth
#'
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @importFrom SummarizedExperiment rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom methods is
#'
#' 
#' 
#' @examples
#' ## DO NOT RUN.  USE THE VIGNETTE CODE TO RUN YOUR ANALYSIS.
#' library(SingleCellExperiment)
#' data(binary_tree_1)
#' sce <- SingleCellExperiment(assays = list(counts=t(binary_tree_1$raw_data),
#' logcounts = t(binary_tree_1$normalized_data)))
#' sce <- PrepareTotem(sce)
#' sce <- RunDimRed(sce,dim.red.method="lmds",
#' dim.reduction.par.list=list(ndim=5))
#' sce <- RunClustering(sce,N.clusterings=200) # Use default N.clusterings=10000
#' sce <- SelectClusterings(sce,selection.method = 1,selection.N.models = 2)
#' sce <- RunSmoothing(sce)
#' clustering_name <-  ReturnTrajNames(sce)[1]
#' dynwrap_object <- ReturnSlingshotObject(sce,clustering_name)

ReturnSlingshotObject.SingleCellExperiment <- function(object,traj.name) {
  
  
  return(metadata(object)$totem$slingshot_trajectory[[traj.name]])
  
  
  
}

#' @rdname ReturnSlingshotObject
#' @aliases ReturnSlingshotObject
setMethod("ReturnSlingshotObject", signature(object = "SingleCellExperiment"),
          ReturnSlingshotObject.SingleCellExperiment)










