#' An example dataset with a tree-shaped trajectory
#'
#' @docType data
#'
#' @usage data(binary_tree_1)
#'
#' @format raw_data, raw_normalized and clustering
#'
#' @keywords datasets
#'
#' @source \url{https://doi.org/10.5281/zenodo.1443566}
#'
#' @examples
#' data(binary_tree_1)
"binary_tree_1"


#' @title Find lineages of an MST from a starting cluster
#'
#' @description
#' The function finds all lineages of an MST. The starting cluster needs
#' to be provided.
#'
#' @param mst.out an object of class "mst"
#' @param start.cluster a character or numeric denoting the starting cluster.

#' @return a list that includes the lineages
#'
#' @keywords lineage MST
#'
#'
FindLineages <- function(mst.out,start.cluster)
{
  
  lineages <- start.cluster
  while (TRUE)
  {
    lineages_new <- lineages
    for (lineage in lineages)
    {
      lineage_split <- unlist(strsplit(lineage,"-"))
      last_cluster <- lineage_split[length(lineage_split)]
      last_cl_cons <- mst.out[last_cluster,]
      adjacent_clusters <- names(last_cl_cons[last_cl_cons != 0])
      for (cluster in adjacent_clusters)
      {
        if (!grepl(paste0("-",cluster,"-"),lineage)&&cluster!=start.cluster)
        {
          lineages_new <- c(lineages_new,paste0(lineage,"-",cluster))
          lineages_new <- lineages_new[lineages_new != lineage]
        }
      }
    }
    if (paste(lineages,collapse = "_") == paste(lineages_new,collapse = "_"))
    {
      break
    }
    lineages <- lineages_new
  }
  return(lineages)
}


#' @title Run k-medoids (CLARA) clustering
#'
#' @description
#' The function runs CLARA algorithm multiple times, with the 
#' number of clusters (k) corresponding to the parameter \code{ks}, 
#' a character vector of unique k values (e.g. \code{c("9", "9.1", "9.2")}).
#' @param X a matrix that specifies the embedding (dim.reduced data).
#' @param ks a character or numeric denoting the starting cluster.
#' @param min.cluster.size a positive integer specifying 
#' the minimum number of cells the clusters can have. Clusterings that don't
#' meet this requirement are filtered out. Note that this should be generally
#' greater or equal to the number of dimensions 
#' @param verbose a logical for denoting whether the index of the clustering
#' that is generated should be printed for progress following purposes.
#'
#' @importFrom cluster clara
#'
#' @return a list that includes the clustering results
#'
#' @keywords clara cluster
#'
#'
RunCLARA <- function(X,ks,min.cluster.size,verbose)
{
  clusterings <- list()
  j <- 1
  for (ks_ in ks) {
    if (verbose) {
      message("clustering: ",j,"/",length(ks))
      j <- j + 1
    }
    clusterings[[ks_]] <- clara(X,k = as.numeric(gsub("\\..*","",ks_)),
                                rngR = TRUE,cluster.only = TRUE)
  }
  names(clusterings) <- ks
  
  for (k in ks) {
    clustering <- clusterings[[k]]
    if (min(table(clustering)) < min.cluster.size) {
      clusterings[[k]] <- NULL
    }
  }
  ks <- names(clusterings)
  
  message("number of clustering results: ", length(ks))
  
  return(clusterings)
}



#' @title Generate MSTs from clustering results
#'
#' @description
#' The function generates MSTs using the Slingshot distance method.
#' @param X a matrix that specifies the embedding (dim.reduced data).
#' @param clusterings a named list of clustering results.
#'
#' @importFrom ape mst
#'
#' @return a list of MSTs
#'
#' @keywords ape MST
#'
#'
GenerateMST <- function(X,clusterings)
{
  ks <- names(clusterings)
  msts <- list()
  for (k in ks) {
    clustering <- clusterings[[k]]
    centroids <- ClusterCentroids(clustering,X)
    mst_out <- mst(SlingshotDistanceMatrix(x=X,
                                           centers=centroids,
                                           clusters=clustering))
    msts[[k]] <- mst_out
  }  
  return(msts)
}



#' @title Calculate connectivity
#'
#' @description
#' The function calculates the node-specific connectivity 
#' for each MST in a list.
#' @param clusterings a named list of clustering results.
#' @param mst.list a named list of MSTs
#'
#' @importFrom plyr mapvalues
#'
#' @return a list of MSTs
#'
#' @keywords ape MST
#'
#'
Connectivity <- function(clusterings,mst.list)
{
  ks <- names(mst.list)
  connectivities <- list()
  for (k in ks) {
    clustering <- clusterings[[k]]
    mst_out <- mst.list[[k]]
    # creates the max-scaled connectivity
    a <- apply(mst_out,1,function(x) sum(x!=0))
    b <- a/(max(apply(mst_out,1,function(x) sum(x!=0))))
    d <- mapvalues(clustering,names(b),b)
    connectivities[[k]] <- d
  }
  
  return(connectivities)
}



#' @title Calculate cluster centroids
#'
#' @description
#' The function calculates the cluster centroids.
#' @param clustering a vector of cluster labels.
#' @param X a matrix that specifies the embedding (dim.reduced data).
#'
#' @return a matrix
#'
#' @keywords cluster centroids center
#'
#'
ClusterCentroids <- function(clustering,X)
{
  centroids <- do.call(rbind,
                       lapply(split(seq_len(nrow(X)),clustering),
                              function(x) apply(X[x,,drop=FALSE],2,mean)))
  
  return(centroids)
}






#' @title Create milestone network from lineages
#'
#' @description
#' The function creates a milestone network from lineages.
#' @param lineages a list of vectors that specify the lineages.
#'
#' @return a data frame
#'
#' @keywords milestone network lineage
#'
#'
CreateMilestoneNetwork <- function(lineages)
{
  milestone_network <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("from", "to", "length","directed")
  colnames(milestone_network) <- x
  for (lineage in lineages)
  {
    milestones <- unlist(strsplit(lineage ,"-"))
    milestone_1 <- 1
    milestone_2 <- 2
    while(milestone_1 != length(milestones))
    {
      # dynverse requires length. here, set to 1 for all milestone pairs
      milestone_network_ <- c(from=milestones[milestone_1],
                              to=milestones[milestone_2],
                              length=1,
                              directed=FALSE)
      milestone_network <- rbind(milestone_network,milestone_network_)
      
      milestone_1 <- milestone_1 + 1
      milestone_2 <- milestone_2 + 1
      
    }
  }
  x <- c("from", "to", "length","directed")
  colnames(milestone_network) <- x
  milestone_pairs <- paste0(milestone_network[,1],"-",milestone_network[,2])
  milestone_network <- milestone_network[!duplicated(milestone_pairs),]
  milestone_network$length <- as.numeric(milestone_network$length)
  milestone_network$directed <- as.logical(milestone_network$directed)
  
  return(milestone_network)
}



#' @title Create milestone progressions
#'
#' @description
#' The function creates the progressions for the dynwrap object.
#' @param clustering a vector of cluster labels.
#' @param X a matrix that specifies the embedding (dim.reduced data).
#' @param milestone.network a data frame of the milestone network.
#'
#' @return a data frame
#'
#' @keywords milestone progressions
#'
#'
CreateProgressions <- function(clustering,X,milestone.network)
{
  progressions <- list()
  clustering_ <- clustering
  names(clustering_) <- rownames(X)
  for (i in seq_len(nrow(X)))
  {
    cell <- rownames(X)[i]
    
    cluster <- clustering_[cell]

    tos <- milestone.network[milestone.network$from==cluster,"to"]
    
    # cell is at an end point cluster
    if (identical(tos, character(0)))
    {
      to <- cluster
      from <- sample(milestone.network[milestone.network$to==to,"from"],1)
      prog <- 1
    } else {
      from <- cluster
      to <- sample(milestone.network[milestone.network$from==from,"to"],1)
      prog <- 0
    }
    
    progression <- c(cell,
                     from,
                     to,
                     prog
    )
    progressions[[cell]] <- progression
  }
  
  progressions <- as.data.frame(do.call(rbind,progressions))
  colnames(progressions) <- c("cell_id","from","to","percentage")
  progressions$percentage <- as.numeric(progressions$percentage)
  rownames(progressions) <- seq_len(nrow(progressions))
  
  return(progressions)
}





#' @title Visualize an MST
#'
#' @description
#' The function visualizes an MST and the underlying clustering.
#' @param viz.dim.red a matrix of the visualization.
#' @param X a matrix that specifies the embedding (dim.reduced data).
#' @param clustering a vector of cluster labels.
#' @param clustering.name a character, the name of the clustering.
#' @param return.traj a logical whether to return the trajectory as well.
#'
#' @importFrom dynplot plot_dimred
#' @importFrom ggplot2 ggtitle theme element_text
#' @importFrom ape mst
#' @importFrom dynwrap wrap_data add_trajectory
#'
#' @return a ggplot2 object
#'
#' @keywords MST visualize
#'
#'
VizOneMST <- function(viz.dim.red,X,clustering,
                      clustering.name,return.traj=FALSE)
{
  centroids <- ClusterCentroids(clustering,X)
  mst_out <- mst(SlingshotDistanceMatrix(x=X,
                                         centers=centroids,
                                         clusters=clustering))
  start_clusters <- rownames(mst_out)[apply(mst_out,1,function(x)sum(x==1))==1]
  start_cluster <- sample(start_clusters,1)
  lineages <- FindLineages(mst_out,start_cluster)
  milestone_ids <- names(table(clustering))
  
  milestone_network <- CreateMilestoneNetwork(lineages)
  progressions <- CreateProgressions(clustering,X,milestone_network)
  
  traj <- wrap_data(cell_ids = rownames(X)) %>%
    add_trajectory(milestone_ids = as.character(milestone_ids),
                   milestone_network = milestone_network,
                   progressions = progressions)
  
  p <- plot_dimred(traj,label_milestones = TRUE,dimred = viz.dim.red) +
    ggtitle(clustering.name) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (return.traj) {
    return(list(p=p,traj=traj))
  } else {
    return(p)
  }
}







#' @title Run Slingshot smoothing for a clustering
#'
#' @description
#' The function runs the Slingshot smoothing for a base clustering.
#'
#' @param rd a matrix that includes the dimensionally reduced data.
#' @param labels a vector that includes the base clustering.
#' @param start.clus an integer or a character 
#' for the name of the start cluster.
#' @param slingshot.par.list a list that includes the parameters of
#' Slingshot.
#' 
#' @import purrr
#' @import dplyr
#' @import dynwrap
#' @import tibble
#' @import slingshot
#' @importFrom princurve project_to_curve
#' @importFrom utils head
#' @importFrom stats dist
#' 
#' @return a list that includes objects of 
#' classes dynwrap and PseudotimeOrdering
#'
#' @keywords slingshot smoothing
#'
#'
RunSlingshotSmoothing <- function(rd,labels,start.clus,slingshot.par.list) {
  sds <- slingshot(
    rd,
    labels,
    shrink = slingshot.par.list$shrink,
    reweight = slingshot.par.list$reweight,
    reassign = slingshot.par.list$reassign,
    thresh = slingshot.par.list$thresh,
    maxit = slingshot.par.list$maxit,
    stretch = slingshot.par.list$stretch,
    smoother = slingshot.par.list$smoother,
    shrink.method = slingshot.par.list$shrink.method,
    dist.method=slingshot.par.list$dist.method,
    start.clus=start.clus
  )
  start_cell <- apply(slingPseudotime(sds), 1, min) %>% 
    sort() %>% head(1) %>% names()
  start.clus <- labels[[start_cell]]
  centroids <- ClusterCentroids(labels,rd)
  ccdm <- as.matrix(dist(centroids))
  from <- to <- NULL
  lineages <- slingLineages(sds)
  lineage_ctrl <- slingParams(sds)
  cluster_network <- lineages %>%
    map_df(~ tibble(from = .[-length(.)], to = .[-1])) %>%
    unique() %>%
    mutate(
      length = ccdm[cbind(from, to)],
      directed = TRUE
    )
  cluster <- slingClusterLabels(sds)
  lin_assign <- apply(slingCurveWeights(sds), 1, which.max)
  progressions <- map_df(seq_along(lineages), function(l) {
    ind <- lin_assign == l
    lin <- lineages[[l]]
    pst.full <- slingPseudotime(sds, na = FALSE)[,l]
    pst <- pst.full[ind]
    means <- sapply(lin, function(clID){
      stats::weighted.mean(pst.full, cluster[,clID])
    })
    non_ends <- means[-c(1,length(means))]
    edgeID.l <- as.numeric(cut(pst, breaks = c(-Inf, non_ends, Inf)))
    from.l <- lineages[[l]][edgeID.l]
    to.l <- lineages[[l]][edgeID.l + 1]
    m.from <- means[from.l]
    m.to <- means[to.l]
    pct <- (pst - m.from) / (m.to - m.from)
    pct[pct < 0] <- 0
    pct[pct > 1] <- 1
    tibble(cell_id=names(which(ind)),from=from.l,to=to.l,percentage=pct)
  })
  output <-
    dynwrap::wrap_data(
      cell_ids = rownames(rd)
    ) %>%
    dynwrap::add_trajectory(
      milestone_network = cluster_network,
      progressions = progressions
    )
  sds@metadata$slingParams$dist.method <- slingshot.par.list$dist.method
  return(list(dynwrap=output,slingshot=sds))
}





#' @title Calculate the Slingshot distance matrix
#'
#' @description
#' The function calculates the Slingshot distance matrix of a clustering.
#' The method is the covariance-based method (Mahalanobis) used by default 
#' in Slingshot. The function has been copied from 
#' the source code of the Slingshot R package.
#' https://github.com/kstreet13/slingshot
#'
#' @param x a matrix for the embedding
#' @param centers a matrix for the cluster centroids
#' @param clusters a vector for the clustering
#' 
#' @importFrom stats cov
#' 
#' @return a matrix of the cluster distances
#'
#' @keywords distance matrix Slingshot Mahalanobis covariance
#'
#'
SlingshotDistanceMatrix <- function(x,centers,clusters)
{
  dist.method <- "slingshot"
  full <- (dist.method=="scaled.full"||
             (dist.method=="slingshot"&&min(table(clusters))>ncol(x)))
  
  nclust <- nrow(centers)
  output_dimnames <- list(rownames(centers), rownames(centers))
  output <- matrix(0, nclust, nclust, dimnames=output_dimnames)
  
  # Computing the covariances (possibly with weights).
  all.cor <- vector("list", nclust)
  names(all.cor) <- rownames(centers)
  
  if (is.matrix(clusters)) {
    # Treating the weights as effective frequencies.
    clusters <- clusters/rowSums(clusters)
    for (i in seq_along(all.cor)) {
      curweight <- clusters[,i]
      out <- t(t(x) - centers[i,]) * sqrt(curweight)
      all.cor[[i]] <- crossprod(out)/(sum(curweight) - 1)
    }
  } else {
    for (i in seq_along(all.cor)) {
      all.cor[[i]] <- cov(x[which(clusters==names(all.cor)[i]),, drop = FALSE])
    }
  }
  
  for (i in seq_len(nclust)) {
    mu1 <- centers[i,]
    clus1 <- rownames(centers)[i]
    s1 <- all.cor[[i]]
    if (!full) {
      s1 <- diag(diag(s1))
    }
    
    for (j in seq_len(i - 1L)) {
      mu2 <- centers[j,]
      clus2 <- rownames(centers)[j]
      s2 <- all.cor[[j]]
      if (!full) {
        s2 <- diag(diag(s2))
      }
      
      diff <- mu1 - mu2
      d <- sqrt(as.numeric(t(diff) %*% solve(s1 + s2) %*% diff))
      output[i,j] <- output[j,i] <- d
    }
  }
  
  return(output)
}





