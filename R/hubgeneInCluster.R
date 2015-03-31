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
#' @seealso \code{\link{cogena}} and \code{\link{geneInCluster}}
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
#' #hubgeneInCluster
#' hubgeneInCluster(cogena_result, "kmeans", "3", "2")
#' 
#' @export
setGeneric("hubgeneInCluster", function(object, method, nClusters, ith) standardGeneric("hubgeneInCluster"))

#' @rdname hubgeneInCluster
#' @exportMethod hubgeneInCluster
#' @aliases hubgeneInCluster
setMethod("hubgeneInCluster", signature(object="cogena"),
          function (object, method=clusterMethods(object), nClusters=nClusters(object), ith){
              #ith is the ith cluster enquerying
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))
              ith <- match.arg(ith, as.character(seq(1:as.numeric(nClusters))))
              ith <- as.numeric(ith)
              
              geneincluster <- geneInCluster (object, method, nClusters, as.character(ith))
              #geneDist <- amap::Dist(mat(object)[geneincluster,], method=object@metric, nbproc=object@ncore)
              geneDist <- object@Distmat
              geneAdjacency <- (1 - as.matrix(geneDist))
              geneAdjacency <- geneAdjacency[]
              names(which.max(colSums(geneAdjacency)))
          })
