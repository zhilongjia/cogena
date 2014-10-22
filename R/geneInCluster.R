setGeneric("geneInCluster", function(object, ...) standardGeneric("geneInCluster"))
#' @export geneInCluster
setMethod("geneInCluster", signature(object="clena"),
          function (object, method=clusterMethods(object), nClusters=nClusters(object), ith){
              #ith is the ith cluster enquerying
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))
              ith <- match.arg(ith, as.character(seq(1:ith)))
              ith <- as.numeric(ith)
              if (method %in% c("hierarchical", "diana", "agnes")) {
                  gene <- names(which (cutree(clusters(object, method), k=nClusters)==ith))
              } else {
                  gene <- names(which(clusters(object, method)[[nClusters]]$cluster==ith))
              }
              return (gene)
          })