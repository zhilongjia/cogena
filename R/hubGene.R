#' Hub genes in a cluster
#' 
#' Show hub genes with p-values in a cluster.
#' 
#' @param geneC a gene names vector
#' @param score_threshold a threshold for the combined scores of the interactions.
#' @param input_directory a directory to store data downloaded from STRING.
#' @return a vector with p-values
#' @export
#' @import igraph
#' @import STRINGdb
#' @rdname hubGene
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
#' geneC3 <- geneInCluster(clen_res, "hierarchical", "3", "3")
#' res <- hubGene(geneC3, input_directory="~/tmp/STRINGdb")
#' }
#' 
setGeneric("hubGene", function(geneC, score_threshold=0, input_directory="") 
    standardGeneric("hubGene"))

#' @rdname hubGene
#' @aliases hubGene,cluster_methods
setMethod("hubGene", signature(geneC="character"),
    function(geneC, score_threshold=0, input_directory="") {
    suppressWarnings(string_db <- STRINGdb::STRINGdb$new(version="10", species=9606, score_threshold=score_threshold, input_directory=input_directory))
    
    example1_mapped <- string_db$map(as.data.frame(geneC), "geneC", removeUnmappedRows = TRUE, quiet=TRUE)
    hits <- example1_mapped$STRING_id
    p <- string_db$get_subnetwork(hits)
    hubvec <- igraph::hub.score(p)$vector
#     hubgenes <- names(which(hubvec==min(hubvec)))
    rownames(example1_mapped) <- example1_mapped$STRING_id
    names(hubvec) <- example1_mapped[names(hubvec),"geneC"]
    hubvec <- sort(hubvec)
    # hubgenes <- example1_mapped$geneC[which(example1_mapped$STRING_id %in% names(hubvec) )]
    return (hubvec)
})

