#' Plot the Protein-protein interaction for genes in a certain cluster.
#' 
#' Show how genes in a certain cluster are interacted based on String database.
#' please note, it doesnot work if there are more than 400 genes in a cluster.
#'
#' @inheritParams enrichment
#' @param ith the ith cluster.
#' @param score_threshold a threshold for the combined scores of the
#'  interactions. Default is 0.
#' @param ... parameter to plot_network from STRINGdb
#' @return a figure
#' @rdname plotPPI
#' @import STRINGdb
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
#' #plotPPI
#' \dontrun{
#' plotPPI(cogena_result, "kmeans", "3", "3")
#' }
#' 
setGeneric("plotPPI", function(object, method, nClusters, ith, score_threshold=0, ...) 
    standardGeneric("plotPPI"))

#' @rdname plotPPI
#' @aliases plotPPI,cogena_methods
setMethod("plotPPI", signature(object="cogena"),
          function (object, method=clusterMethods(object), 
                    nClusters=nClusters(object), ith,
                    score_threshold=0, ...){
    geneC <- geneInCluster(object, method, nClusters, as.character(ith))
    suppressWarnings(string_db <- STRINGdb$new(version="10", species=9606, score_threshold=score_threshold))
    
    example1_mapped <- string_db$map(as.data.frame(geneC), "geneC", removeUnmappedRows = TRUE, quiet=TRUE)
    hits <- example1_mapped$STRING_id
    string_db$plot_network(hits, ...)
    }
)

