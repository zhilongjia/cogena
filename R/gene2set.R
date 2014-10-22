#' generate relationship between genes and gene sets
#' 
#' Generate relationship between genes (SYMBOL) and gene sets, such as Pathway or GO.
#' 
#' The gene sets could be gmt format. Options could be c2.cgp.v4.0.symbols.gmt,
#' c2.cp.reactome.v4.0.symbols.gmt,  c5.bp.v4.0.symbols.gmt, 
#' c2.cp.kegg.v4.0.symbols.gmt  c2.cp.v4.0.symbols.gmt, c6.all.v4.0.symbols.gmt,
#' msigdb.v4.0.symbols.gmt. Or a user-defined gene sets gmt file.
#' 
#' @param annoList gene set annotation gmt file. The default is c5.bp.v4.0.symbols.gmt from MSigDB Collections.
#' @param annotation an annotation list.
#' @param genenames a SYMBOL gene names charactic vector.
#' @param TermFreq a threshold for the Term Frequence.
#' @return an annotaion matrix.
#' @examples ###
#' @export
#' @rdname gene2set
gene2set <- function(anno, genenames, TermFreq=0) {
    if (is.null(anno)) {
        anno <- "c2.cp.v4.0.symbols.gmt"}
    anno <- match.arg(anno, c("c2.cgp.v4.0.symbols.gmt", "c2.cp.kegg.v4.0.symbols.gmt", 
        "c2.cp.reactome.v4.0.symbols.gmt", "c2.cp.v4.0.symbols.gmt", "c5.bp.v4.0.symbols.gmt", 
        "c6.all.v4.0.symbols.gmt"))
    data(annotation, package="clena")
    annoMarix = annotationListToMatrix(get(anno), genenames)
    anno = subset(annoMarix, select=names(which( colSums(annoMarix)/sapply(anno, length)>=TermFreq )))
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

