#' enrichment for genes in a certain cluster.
#' 
#' A p-value representing the enrichment in interactions of the list of
#' proteins (i.e. the probability to obtain such a number of interactions
#' by chance). It is mainly based on get_ppi_enrichment from STRINGdb.
#'
#' @inheritParams geneInCluster
#' @param score_threshold a threshold for the combined scores of the
#'  interactions. Default is 0.
#' @return p-value
#' @rdname PPIenrichment
#' @import STRINGdb
#' @export
#' @seealso \code{\link{clEnrich}}, \code{\link[STRINGdb]{get_ppi_enrichment}}
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
#' #PPIenrichment
#' PPIenrichment(clen_res, "kmeans", "3", "3")
#' }
#' 
setGeneric("PPIenrichment", function(object, method, nCluster, ith, score_threshold=0) 
    standardGeneric("PPIenrichment"))

#' @rdname PPIenrichment
#' @aliases PPIenrichment,cogena_methods
setMethod("PPIenrichment", signature(object="cogena"),
          function (object, method=clusterMethods(object), 
                    nCluster=nClusters(object), ith,
                    score_threshold=0){
    geneC <- geneInCluster(object, method, nCluster, as.character(ith))
    suppressWarnings(string_db <- STRINGdb$new(version="10", species=9606, score_threshold=score_threshold))
    
    example1_mapped <- string_db$map(as.data.frame(geneC), "geneC", removeUnmappedRows = TRUE, quiet=TRUE)
    hits <- example1_mapped$STRING_id
    string_db$get_ppi_enrichment(hits)$enrichment
    }
)
