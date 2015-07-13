#' Get gene names in a certain cluster.
#' 
#' Get gene names in a certain cluster. This is helpful if user want to get the
#' detail of a cluster.
#'
#' @param object a cogena object
#' @param method a clustering method
#' @param nCluster cluster number
#' @param ith the i-th cluster (should no more than nCluster)
#' @return a character vector containing the gene names.
#' @rdname geneInCluster
#' @export
#' @seealso \code{\link{clEnrich}}
#' @examples
#' data(PD)
#' annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt", 
#' package="cogena")
#' 
#' \dontrun{
#' genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
#' 
#' #summay this cogena object
#' summary(clen_res)
#' 
#' #geneInCluster
#' g1 <- geneInCluster(clen_res, "kmeans", "3", "2")
#' 
#' #Up or Down genes with setting nCluster as "2".
#' g2 <- geneInCluster(clen_res, "kmeans", "2", "1")
#' }
#' 
setGeneric("geneInCluster", function(object, method, nCluster, ith) 
    standardGeneric("geneInCluster"))

#' @rdname geneInCluster
#' @aliases geneInCluster,cogena_methods
setMethod("geneInCluster", signature(object="cogena"),
    function (object, method=clusterMethods(object), 
        nCluster=nClusters(object), ith){
    #ith is the ith cluster enquerying
    method <- match.arg(method, clusterMethods(object))
    nCluster <- match.arg(nCluster, as.character(nClusters(object)))

    ith <- match.arg(ith, as.character(seq(1: as.numeric(nCluster))))
    ith <- as.numeric(ith)

    cluster_size <- geneclusters(object, method, nCluster)
    names(cluster_size) <- rownames(object@mat)
    names(which(cluster_size==ith))
    }
)

