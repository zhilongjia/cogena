#' Get gene names in a certain cluster.
#' 
#' Get gene names in a certain cluster. This is helpful if user want to get the
#' detail of a cluster.
#'
#' @inheritParams enrichment
#' @param ith the ith cluster.
#' @return a character vector containing the gene names.
#' @rdname geneInCluster
#' @export
#' @seealso \code{\link{cogena}}
#' @examples
#' \dontrun{
#' #summay this cogena object
#' summary(cogena_result)
#' 
#' #geneInCluster
#' geneInCluster(cogena_result, "kmeans", "3", "2")
#' }
setGeneric("geneInCluster", function(object, method, nClusters, ith) standardGeneric("geneInCluster"))

#' @rdname geneInCluster
#' @aliases geneInCluster,cogena_methods
setMethod("geneInCluster", signature(object="cogena"),
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
