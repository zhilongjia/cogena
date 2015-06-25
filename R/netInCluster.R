#' Show network in a certain cluster.
#' 
#' It's based on igraph package.
#'
#' @param object a cogena object
#' @param cutoff similarity cut-off. Default is 0.7
#' @param gn A charactor vector from geneInCluster
#' @param ... other parameter to igraph::plot.graph
#' @return a figure.
#' @rdname netInCluster
#' @import igraph
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
#' #geneInCluster
#' g1 <- geneInCluster(cogena_result, "kmeans", "3", "3")
#' 
#' netInCluster(cogena_result, g1)
#' 
setGeneric("netInCluster", function(object, gn, cutoff=0.7, ...) 
    standardGeneric("netInCluster"))

#' @rdname netInCluster
#' @aliases netInCluster,cogena_methods
setMethod("netInCluster", signature(object="cogena"),
          function (object, gn, cutoff=0.7, ...){
              distmat <- 1- as.matrix(object@Distmat)[gn,gn]
              distmat[distmat<=cutoff] <- 0
              net <- igraph::graph.adjacency(distmat, mode="undirected", weighted=TRUE)
              net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)
              plot(net, ...)
          }
)
