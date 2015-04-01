#' Get gene names in each clusters and the expression profiling.
#' 
#' Get gene names in each clusters and the expression profiling. This output 
#' is helpful if user want to analyse the data for other application.
#'
#' @inheritParams enrichment
#' @return a list containing a matrix of cluster_id with expression profiling 
#' and label a vector of the sample labels.
#' @rdname geneExpInCluster
#' @export
#' @seealso \code{\link{cogena}}
#' @examples 
#' data(PD)
#' annofile <- system.file("extdata", "c2.cp.kegg.v4.0.symbols.gmt", 
#' package="cogena")
#' cogena_result <- cogena(DEexprs, nClust=2:3, 
#' clMethods=c("hierarchical","kmeans"), metric="correlation", 
#' method="complete",  annofile=annofile, sampleLabel=sampleLabel, 
#' ncore=1, verbose=TRUE)
#' #summay this cogena object
#' summary(cogena_result)
#' 
#' #geneExpInCluster
#' geneExpInCluster(cogena_result, "kmeans", "3")
#' 
#' 
setGeneric("geneExpInCluster", function(object, method, nClusters) 
    standardGeneric("geneExpInCluster"))

#' @rdname geneExpInCluster
#' @aliases geneExpInCluster,cogena_methods
setMethod("geneExpInCluster", signature(object="cogena"),
    function (object, method=clusterMethods(object), 
        nClusters=nClusters(object)){

    method <- match.arg(method, clusterMethods(object))
    nClusters <- match.arg(nClusters, as.character(nClusters(object)))


    if (method %in% c("hierarchical", "diana", "agnes")) {
        cluster_id <- cutree(clusters(object, method), k=nClusters)
    } else {
        cluster_id <- clusters(object, method)[[nClusters]]$cluster
    }

    clusterGeneExp <- cbind(cluster_id, object@mat)
    clusterGeneExp <- clusterGeneExp[order(clusterGeneExp[,"cluster_id"]),]

    label <- object@sampleLabel
    list(clusterGeneExp=clusterGeneExp, label=label)
    }
)

