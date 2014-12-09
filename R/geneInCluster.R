#' Show gene names in certain cluster.
#' 
#' Show gene names in certain cluster.
#'
#' @inheritParams enrichment
#' @param ith the ith cluster.
#' @return a character vector.
#' @docType methods
#' @rdname geneInCluster
#' 
#' @seealso \code{\link{clena}}
#' @examples 
#' #summay this clena object
#' summary(GSE7621.SAM.clena.cluster)
#' 
#' #geneInCluster
#' geneInCluster(GSE7621.SAM.clena.cluster, "sota", "12", "2")
#' 
#' @export
setGeneric("geneInCluster", function(object, method, nClusters, ith, ...) standardGeneric("geneInCluster"))


#' @exportMethod
#' @aliases geneInCluster
setMethod("geneInCluster", signature(object="clena"),
          function (object, method=clusterMethods(object), nClusters=nClusters(object), ith){
              #ith is the ith cluster enquerying
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))

              ith <- match.arg(ith, as.character(seq(1: as.numeric(nClusters))))
              ith <- as.numeric(ith)
              if (method %in% c("hierarchical", "diana", "agnes")) {
                  cluster_size <- cutree(clusters(object, method), k=nClusters)
              } else {
                  cluster_size <- clusters(object, method)[[nClusters]]$cluster
              }
              names(cluster_size) <- rownames(object@mat)
              names(which(cluster_size==ith))
          })