#' Show hub gene names in certain cluster.
#' 
#' Show hub gene names in certain cluster.
#'
#' @inheritParams enrichment
#' @param ith the ith cluster.
#' @return a character vector.
#' @docType methods
#' @rdname hubgeneInCluster
#' 
#' @seealso \code{\link{clEnrich}} and \code{\link{geneInCluster}}
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
#' #hubgeneInCluster
#' hubgeneInCluster(clen_res, "kmeans", "3", "2")
#' }
#'
#' @export
setGeneric("hubgeneInCluster", 
    function(object, method, nCluster, ith) 
    standardGeneric("hubgeneInCluster"))

#' @rdname hubgeneInCluster
#' @exportMethod hubgeneInCluster
#' @aliases hubgeneInCluster
setMethod("hubgeneInCluster", signature(object="cogena"),
    function (object, method=clusterMethods(object),
        nCluster=nClusters(object), ith){
    #ith is the ith cluster enquerying
    method <- match.arg(method, clusterMethods(object))
    nCluster <- match.arg(nCluster, as.character(nClusters(object)))
    ith <- match.arg(ith, as.character(seq(1:as.numeric(nCluster))))
    ith <- as.numeric(ith)
    
    geneincluster <- geneInCluster(object, method, nCluster, as.character(ith))
    geneDist <- object@Distmat
    geneAdjacency <- (1 - as.matrix(geneDist))
    geneAdjacency <- geneAdjacency[]
    names(which.max(colSums(geneAdjacency)))
    }
)

