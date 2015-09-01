#' generate relationship between genes and gene-sets
#' 
#' Generate relationship between genes (gene SYMBOL) and gene-sets, 
#' such as Pathway or GO.
#' 
#' @inheritParams gmt2list
#' @param genenames a SYMBOL gene names charactic vector.
#' @param TermFreq a threshold for the Term Frequence. Default is zero.
#' @param annotation a value returned by \code{\link{gmt2list}}.
#' @return an gene and gene-set relationship matrix 
#' @examples
#' data(Psoriasis)
#' 
#' #annotaion
#' annoGMT <- "c2.cp.kegg.v5.0.symbols.gmt.xz"
#' annofile <- system.file("extdata", annoGMT, package="cogena")
#' # the DEG gene-sets matrix
#' anno <- gene2set(annofile, rownames(DEexprs))
#' @export
#' @rdname gene2set
gene2set <- function(annofile=NULL, genenames, TermFreq=0) {
    if (is.null(annofile)) {
        annofile <- system.file("extdata", "c2.cp.v5.0.symbols.gmt.xz", 
            package="cogena")
    }

    annoList <- gmt2list(annofile)
    annoMatrix <- annotationListToMatrix(annoList, genenames)
    anno=annoMatrix[,
        names(which( colSums(annoMatrix)/ncol(annoMatrix)>=TermFreq ))]
    return(anno)
}


#' @rdname gene2set
annotationListToMatrix <- function(annotation, genenames) {
    annotation.matrix <- matrix(FALSE, ncol=length(annotation), 
        nrow=length(genenames))
    colnames(annotation.matrix) <- names(annotation)
    rownames(annotation.matrix) <- genenames
    for ( i in 1:length(annotation) ){
        annot <- names(annotation)[i]
        genes <- as.character(annotation[[i]][!is.na(annotation[[i]])])
        genes.common <- intersect(genenames, genes)
        annotation.matrix[genes.common, annot ] <- TRUE
    }
    return(annotation.matrix)
}

