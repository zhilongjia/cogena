#' Significance of Gene sets enrichment.
#' 
#' Caculating the significance of Gene sets enrichment based on the 
#' hypergeometric test. This function is mainly used internally.
#' 
#' Here the genes in annotation can be a varity of types. like all the DEG, 
#' up-regualted genes or genes in a cluster. the gene names should be 
#' consistent with the genes in the gene sets.
#' 
#' 
#' @param genenames a vector of gene names.
#' @param annotation data.frame with the gene (like all the differentially 
#' expressed genes) in row, gene set in column.
#' @param annotationGenesPop data.frame with the gene in row, gene set 
#' in column. Here genes are genes in population with filering the 
#' non-nformative genes better.
#' @return a vector with P-values.
#' @export
#' @examples
#' data(Psoriasis)
#' data(AllGeneSymbols)
#' annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt.xz", package="cogena")
#' annoBG <- gene2set(annofile, AllGeneSymbols)
#' res <- PEI(rownames(DEexprs)[1:200], gene2set(annofile, rownames(DEexprs)[1:200]), annoBG)
#' 

PEI <- function(genenames, annotation, annotationGenesPop) {

    # presteps for calculating p-value
    #pei <- vector("numeric", ncol(annotation))
    pei <- vector(mode="numeric")
    NumGenes <- nrow(annotationGenesPop)
    for (j in colnames(annotation) ){
        NumGenesWithinPathway <- sum(annotationGenesPop[,j])
        NumGenesWithinCluster <- length(genenames)
        NumGenesWithinCluster.Pathway <- 
            sum(annotation[genenames,j]) -1 #minus one due to phyper
        if (NumGenesWithinPathway<2 | 
            NumGenesWithinCluster<2 | 
            NumGenesWithinCluster.Pathway<2){
            pei[j] = NA
        } else {
            pei[j] <- phyper(NumGenesWithinCluster.Pathway, 
                NumGenesWithinPathway, NumGenes-NumGenesWithinPathway, 
                NumGenesWithinCluster, lower.tail=FALSE)
        }
    }
    pei
}
