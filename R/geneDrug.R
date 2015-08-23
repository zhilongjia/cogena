#' The drug and gene interaction
#' 
#' Get the drug and gene interaction from DGIdb.
#' 
#' @param genenames the gene names. A vector.
#' @return a data.frame showing drugs and the interactionTypes
#' @export
#' @rdname geneDrug
#' @import jsonlite
#' @examples 
#' 
#' data(Psoriasis)
#' annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt.xz", 
#' package="cogena")
#' 
#' \dontrun{
#' genecl_result <- coExp(DEexprs, nClust=c(2,10), clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
#' 
#' geneC3 <- geneInCluster(clen_res, "hierarchical", "10", "3")
#' 
#' dgi <- drugGene(geneC3)
#' }
#' 
#' 
setGeneric("geneDrug", function(genenames) 
    standardGeneric("geneDrug"))

#' @rdname geneDrug
#' @aliases geneDrug
setMethod("geneDrug", signature(genenames="character"), 
    function(genenames) {
        genenames <- paste0(genenames, collapse=",")
        dgiAPI <- "http://dgidb.genome.wustl.edu/api/v1/interactions.json?genes="
        dgi_query <- paste0(dgiAPI, genenames)
    
        res <- jsonlite::fromJSON(dgi_query)$matchedTerms
    
        return (res)
})

