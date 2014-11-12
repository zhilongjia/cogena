#' generate relationship between genes and gene sets
#' 
#' Generate relationship between genes (SYMBOL) and gene sets, such as Pathway or GO.
#' 
#' The gene sets could be gmt format. Options could be c2.cgp.v4.0.symbols.gmt,
#' c2.cp.reactome.v4.0.symbols.gmt,  c5.bp.v4.0.symbols.gmt, 
#' c2.cp.kegg.v4.0.symbols.gmt  c2.cp.v4.0.symbols.gmt, c6.all.v4.0.symbols.gmt,
#' msigdb.v4.0.symbols.gmt. Or a user-defined gene sets gmt file.
#' 
#' @param anno a string or a gmt file. The default is "c2.cp.v4.0.symbols.gmt"
#' from MSigDB Collections. Others like "c2.cgp.v4.0.symbols.gmt", 
#' "c2.cp.kegg.v4.0.symbols.gmt", "c2.cp.reactome.v4.0.symbols.gmt", 
#' "c5.bp.v4.0.symbols.gmt", "c6.all.v4.0.symbols.gmt". Or a gmt format file.
#' @param genenames a SYMBOL gene names charactic vector.
#' @param TermFreq a threshold for the Term Frequence. Default is zero.
#' 
#' @examples
#' data(PD)
#' annoGMT <- "c2.cp.kegg.v4.0.symbols.gmt"
#' annoGMT <- system.file("data", annoGMT, package="clena")
#' anno = gene2set(anno=annoGMT, difSAM_GSE7621_DEG, TermFreq=0)
#' annotationGenesPop = gene2set(anno=annoGMT, rownames(GSE7621.filtered.expr), TermFreq=0)
#' annotationGenesPop <- annotationGenesPop[,colnames(anno)]
#' 
#' @export
#' @rdname gene2set
gene2set <- function(anno, genenames, TermFreq=0) {
    if (is.null(anno)) {
        anno <- "c2.cp.v4.0.symbols.gmt"
        annofile <- system.file("data", annoGMT, package="clena")
    }
    
    annoList <- GSEABase::geneIds(GSEABase::getGmt(anno))
    annoMatrix <- annotationListToMatrix(annoList, genenames)
    # annoMatrix <- annoMatrix[which(apply(annoMatrix, 1, sum)>0),]
    # if (nrwo(annoMatrix) < length(genenames)){
    #     warning("Some genes does NOT exist in the gene sets. Please revise the expression profilings.") 
    # }
    anno = subset(annoMatrix, select=names(which( colSums(annoMatrix)/sapply(anno, length)>=TermFreq )))
}



annotationListToMatrix <- function(annotation, genenames) {
    annotation.matrix <- matrix(FALSE, ncol=length(annotation), nrow=length(genenames))
    colnames(annotation.matrix) <- names(annotation)
    rownames(annotation.matrix) <- genenames
    for ( i in 1:length(annotation) )
    {
        annot <- names(annotation)[i]
        genes <- as.character(annotation[[i]][!is.na(annotation[[i]])])
        genes.common <- intersect(genenames, genes)
        annotation.matrix[genes.common, annot ] <- TRUE
    }
    return(annotation.matrix)
}

