#' generate relationship between genes and gene-sets
#' 
#' Generate relationship between genes (gene SYMBOL) and gene-sets, such as Pathway or GO.
#' 
#' @inheritParams gmt2list
#' @param genenames a SYMBOL gene names charactic vector.
#' @param TermFreq a threshold for the Term Frequence. Default is zero.
#' 
#' @examples
#' please ?cogena
#' @export
#' @rdname gene2set
gene2set <- function(anno, genenames, TermFreq=0) {
    if (is.null(anno)) {
        anno <- "c2.cp.v4.0.symbols.gmt"
        annofile <- system.file("data", annoGMT, package="cogena")
    }

    annoList <- gmt2list(anno)
    annoMatrix <- annotationListToMatrix(annoList, genenames)
    anno = subset(annoMatrix, select=names(which( colSums(annoMatrix)/sapply(anno, length)>=TermFreq )))
}


#' @rdname gene2set
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

