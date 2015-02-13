#' Correlation in the cluster of a cogena object
#' 
#' Correlation in the cluster of a cogena object. This is helpful if the number 
#' of genes in cluster are small.
#' @inheritParams clusterMethods
#' @inheritParams enrichment
#' @inheritParams geneInCluster
#' @param corMethod a character string indicating which correlation coefficient 
#' (or covariance) is to be computed. One of "pearson" (default), "kendall", or 
#' "spearman", can be abbreviated.
#' @param plotMethod Character, the visualization method of correlation matrix 
#' to be used. Currently, it supports seven methods, named "circle" (default), 
#' "square", "ellipse", "number", "pie", "shade" and "color". See examples in 
#' \code{\link[corrplot]{corrplot}} for details
#' @param type Character, "full" (default), "upper" or "lower", display full 
#' matrix, lower triangular or upper triangular matrix. See examples in 
#' \code{\link[corrplot]{corrplot}} for details
#' @param ... other parameters to \code{\link[corrplot]{corrplot}} function.
#' @rdname corInCluster
#' @seealso \code{\link{cogena}} \code{\link[corrplot]{corrplot}}
#' @export
#' @examples
#' data(PD)
#' corInCluster(cogena_result, "kmeans", "8", "8")
#' corInCluster(cogena_result, "kmeans", "8", "8", plotMethod="square")

setGeneric("corInCluster", function(object, method, nClusters, ith, corMethod="pearson", ...) standardGeneric("corInCluster"))

#' @aliases corInCluster,cogena_methods
setMethod("corInCluster", signature(object="cogena"), 
          function (object, method=clusterMethods(object), 
                    nClusters=nClusters(object), ith,
                    corMethod="pearson", plotMethod = "circle", type = "upper",...){
              
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))
              geneExp <- geneExpInCluster(object, method, nClusters)$clusterGeneExp
              geneExpCluster <- geneExp[geneExp[,"cluster_id"] == ith,-1]
              M <- cor(t(geneExpCluster), method=corMethod)
              corrplot::corrplot(M, method = plotMethod, type = type, ...)
    
})