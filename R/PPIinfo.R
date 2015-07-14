#' Protein-protein interaction infromation for genes
#' 
#' 
#' PPI information for genes in a certain cluster. It will return a short link
#'  to the network page of our STRING website that shows the protein 
#'  interactions between input genes and a p-value representing the enrichment 
#'  in interactions of the list of vproteins (i.e. the probability to obtain
#'   such a number of interactions by chance). It is mainly based on 
#'   get_link and get_summary from STRINGdb.
#'
#' @inheritParams hubGene
#' @param ... to get_link from STRINGdb
#' @return summary of the network constructed by input genes.
#' @rdname PPIinfo
#' @import STRINGdb
#' @export
#' @seealso \code{\link{clEnrich}}, \code{\link[STRINGdb]{get_link}}
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
#' #PPIinfo
#' geneC3 <- geneInCluster(clen_res, "hierarchical", "3", "3")
#' PPIinfo(geneC3, input_directory="~/tmp/STRINGdb")
#' }
#' 
setGeneric("PPIinfo", function(geneC, score_threshold=0, input_directory="", ...) 
    standardGeneric("PPIinfo"))

#' @rdname PPIinfo
#' @aliases PPIinfo,cluster_methods
setMethod("PPIinfo", signature(geneC="character"),
          function(geneC, score_threshold=0, input_directory="", ...){

    suppressWarnings(string_db <- STRINGdb$new(version="10", species=9606, score_threshold=score_threshold, input_directory=input_directory))
    
    example1_mapped <- string_db$map(as.data.frame(geneC), "geneC", removeUnmappedRows = TRUE, quiet=TRUE)
    hits <- example1_mapped$STRING_id
    if (length(hits) <= 400){
        url <- string_db$get_link(hits, ...)
        net_summary <- string_db$get_summary(hits)
        cat ("STRING Website:", url, "\n",  net_summary)
    } else {
        net_summary <- string_db$get_summary(hits)
        cat (net_summary)
    }

})
