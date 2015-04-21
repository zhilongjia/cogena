#' This function creates a SignedGeneSet object from user-specified identifiers.
#'
#' create_SignedGeneSet is revised from gCMAPWeb::create_SignedGeneSet 
#' funcntion. The License is Artistic-2.0.
#'
#' @title SignedGeneSet creator
#' @param genes.up Charactor vector, up-regulated genes.
#' @param genes.down Charactor vector, down-regulated genes.
#' @return SignedGeneSet object
#' @author Thomas Sandmann, Zhilong Jia
#' @importClassesFrom gCMAP SignedGeneSet
#' @importMethodsFrom gCMAP geneSign SignedGeneSet
#' @importMethodsFrom GSEABase SymbolIdentifier
#' @examples 
#' genes.up <- c("ZNF474", "CCDC100", "ANKRD43")
#' genes.down <- c("RNUXA", "TNFAIP8", "MARCH3", "FLJ90650", "PTMAP2")
#' gs <- create_SignedGeneSet(genes.up, genes.down)
#'
create_SignedGeneSet <- function(genes.up=NULL, genes.down=NULL){

    ## remove genes present in both categories
    if( length( genes.up ) > 0 & length( genes.down) >0 ){
        up.and.down <- intersect(genes.up,genes.down)
        if( length( up.and.down) > 0){
            genes.up <- genes.up[-match( up.and.down, genes.up)]
            genes.down <- genes.down[-match( up.and.down, genes.down)]
        }
    }
  
    ## create 'sign' vector
    gene.ids <- c(genes.up, genes.down)
    
    signs <- c(rep("up", length(genes.up)),
             rep("down", length(genes.down)))
    names(signs) <- gene.ids

    ## create SignedGeneSet object
    gs <- SignedGeneSet(gene.ids, setName="queried.genes", geneSign=signs, 
                        geneIdType=GSEABase::SymbolIdentifier())
  return(gs)
}
