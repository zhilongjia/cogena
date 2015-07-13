#' Correlation in the cluster of a cogena object
#' 
#' Correlation in the cluster of a cogena object. This is helpful if the number 
#' of genes in cluster are small.
#' @inheritParams geneInCluster
#' @param corMethod a character string indicating which correlation coefficient 
#' (or covariance) is to be computed. One of "pearson" (default), "kendall", or 
#' "spearman", can be abbreviated.
#' @param plotMethod the visualization method of correlation matrix 
#' to be used. Currently, it supports seven methods, named "circle" (default), 
#' "square", "ellipse", "number", "pie", "shade" and "color". See examples in 
#' \code{\link[corrplot]{corrplot}} for details
#' @param type "full" (default), "upper" or "lower", display full 
#' matrix, lower triangular or upper triangular matrix. See examples in 
#' \code{\link[corrplot]{corrplot}} for details
#' @param ... other parameters to \code{\link[corrplot]{corrplot}} function.
#' @return a correlation figure.
#' @rdname corInCluster
#' @importFrom corrplot corrplot
#' @seealso \code{\link{clEnrich}} \code{\link[corrplot]{corrplot}}
#' @export
#' @examples
#' data(PD)
#' annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt", 
#' package="cogena")
#' genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
#' 
#' corInCluster(clen_res, "kmeans", "3", "3")
#' corInCluster(clen_res, "kmeans", "3", "3", plotMethod="square")
#' 
#' 

setGeneric("corInCluster", 
    function(object, method, nCluster, ith, 
        corMethod="pearson", plotMethod = "circle", type = "upper", ...) 
        standardGeneric("corInCluster"))

#' @rdname corInCluster
#' @aliases corInCluster,cogena_methods
setMethod("corInCluster", signature(object="cogena"), 
    function (object, method=clusterMethods(object), 
        nCluster=nClusters(object), ith,
        corMethod="pearson", plotMethod = "circle", type="upper",...){

    method <- match.arg(method, clusterMethods(object))
    nCluster <- match.arg(nCluster, as.character(nClusters(object)))
    geneExp <- geneExpInCluster(object, method, nCluster)$clusterGeneExp
    geneExpCluster <- geneExp[geneExp[,"cluster_id"] == ith,-1]
    M <- cor(t(geneExpCluster), method=corMethod)
    corrplot(M, method = plotMethod, type = type, ...)
    }
)
