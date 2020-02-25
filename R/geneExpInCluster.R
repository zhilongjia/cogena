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
#' @seealso \code{\link{clEnrich}}
#' @examples 
#' data(Psoriasis)
#' annofile <- system.file("extdata", "c2.cp.kegg.v7.01.symbols.gmt.xz", 
#' package="cogena")
#' 
#' \dontrun{
#' genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
#' 
#' #geneExpInCluster
#' geneExp <- geneExpInCluster(clen_res, "kmeans", "3")
#' }
#' 
setGeneric("geneExpInCluster", function(object, method, nCluster) 
    standardGeneric("geneExpInCluster"))

#' @rdname geneExpInCluster
#' @aliases geneExpInCluster,cogena_methods
setMethod("geneExpInCluster", signature(object="cogena"),
    function (object, method=clusterMethods(object), 
        nCluster=nClusters(object)){

    method <- match.arg(method, clusterMethods(object))
    nCluster <- match.arg(nCluster, as.character(nClusters(object)))


    cluster_id <- geneclusters(object, method,nCluster)
    

    clusterGeneExp <- cbind(cluster_id, object@mat)
    clusterGeneExp <- clusterGeneExp[order(clusterGeneExp[,"cluster_id"]),]

    label <- object@sampleLabel
    list(clusterGeneExp=clusterGeneExp, label=label)
    }
)

