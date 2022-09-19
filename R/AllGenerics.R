

#' @export
setGeneric("RunTotem",signature = "object",
           function(object,
                    dim.red.method="lmds",
                    dim.red.features=NULL,
                    dim.reduction.par.list=list(ndim=5),
                    k.range=3:20,
                    min.cluster.size=5,
                    N.clusterings=10000,
                    selection.method=1,
                    selection.N.models=1,
                    selection.k.range=NULL,
                    selection.stratified=FALSE,
                    prior.clustering=NULL,
                    slingshot.par.list=list(shrink = 1,
                                            reweight = TRUE,
                                            reassign = TRUE,
                                            thresh = 0.001,
                                            maxit = 10,
                                            stretch = 2,
                                            smoother = "smooth.spline",
                                            shrink.method = "cosine",
                                            dist.method="slingshot")) {
             standardGeneric("RunTotem")
           })


#' @export
setGeneric("PrepareTotem",signature = "object",
           function(object) {
             standardGeneric("PrepareTotem")
           })


#' @export
setGeneric("RunDimRed",signature = "object",
           function(object,dim.red.method="lmds",dim.red.features=NULL,
                    dim.reduction.par.list=list(ndim=5)) {
             standardGeneric("RunDimRed")
           })


#' @export
setGeneric("RunClustering",signature = "object",
           function(object,k.range=3:20,min.cluster.size=5,
                    N.clusterings=10000,verbose=FALSE) {
             standardGeneric("RunClustering")
           })


#' @export
setGeneric("SelectClusterings",signature = "object",
           function(object,selection.method=1,selection.N.models=1,
                    selection.k.range=NULL,selection.stratified=FALSE,
                    prior.clustering=NULL) {
             standardGeneric("SelectClusterings")
           })


#' @export
setGeneric("RunSmoothing",signature = "object",
           function(object,
                    slingshot.par.list=list(shrink = 1,
                                            reweight = TRUE,
                                            reassign = TRUE,
                                            thresh = 0.001,
                                            maxit = 10,
                                            stretch = 2,
                                            smoother = "smooth.spline",
                                            shrink.method = "cosine",
                                            dist.method="slingshot")) {
             standardGeneric("RunSmoothing")
           })


#' @export
setGeneric("VizMST",signature = "object",
           function(object,clustering.names=NULL,
                    dim.red.type=NULL,viz.dim.red=NULL) {
             standardGeneric("VizMST")
           })



#' @export
setGeneric("VizClustering",signature = "object",
           function(object,clustering=NULL,dim.red.type=NULL,
                    viz.dim.red=NULL,plot.mst=TRUE) {
             standardGeneric("VizClustering")
           })


#' @export
setGeneric("VizSmoothedTraj",signature = "object",
           function(object,traj.names=NULL,dim.red.type=NULL,
                    viz.dim.red=NULL,plot.pseudotime=FALSE) {
             standardGeneric("VizSmoothedTraj")
           })


#' @export
setGeneric("VizFeatureExpression",signature = "object",
           function(object,traj.name=NULL,feature.names=NULL,dim.red.type=NULL,
                    viz.dim.red=NULL,plot.traj=TRUE) {
             standardGeneric("VizFeatureExpression")
           })


#' @export
setGeneric("VizCellConnectivity",signature = "object",
           function(object,dim.red.type=NULL,viz.dim.red=NULL) {
             standardGeneric("VizCellConnectivity")
           })


#' @export
setGeneric("ReturnTrajNames",signature = "object",
           function(object) {
             standardGeneric("ReturnTrajNames")
           })



#' @export
setGeneric("ReturnClustering",signature = "object",
           function(object,clustering.name) {
             standardGeneric("ReturnClustering")
           })




#' @export
setGeneric("ReturnMSTNetwork",signature = "object",
           function(object,clustering.name) {
             standardGeneric("ReturnMSTNetwork")
           })





#' @export
setGeneric("ReturnSmoothedTrajNetwork",signature = "object",
           function(object,clustering.name) {
             standardGeneric("ReturnSmoothedTrajNetwork")
           })



#' @export
setGeneric("ChangeTrajRoot",signature = "object",
           function(object,traj.name,root.cluster) {
             standardGeneric("ChangeTrajRoot")
           })

#' @export
setGeneric("MergeClusters",signature = "object",
           function(object,clustering.name,merged.clusters,merged.name) {
             standardGeneric("MergeClusters")
           })



#' @export
setGeneric("ReturnDynwrapObject",signature = "object",
           function(object,traj.name) {
             standardGeneric("ReturnDynwrapObject")
           })



#' @export
setGeneric("ReturnSlingshotObject",signature = "object",
           function(object,traj.name) {
             standardGeneric("ReturnSlingshotObject")
           })

