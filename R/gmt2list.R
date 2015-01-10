#' read gmt file as a list
#' 
#' read Gene Matrix Transposed (gmt) file and output a list with the the first column as 
#' the names of items in the list. see \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{Gene Matrix Transposed file format}
#' for more details.
#' 
#' @param anno a gmt file. Examples are from MSigDB Collections.
#' A list of gene set could be find in the vignette of clena.
#'
#' @export
#' 
#' @examples
#' anno <- "c2.cp.kegg.v4.0.symbols.gmt"
#' anno <- system.file("data", anno, package="clena")
#' gmt2list(anno)
#' 
gmt2list <- function(anno){
    x <- scan(anno, what="", sep="\n")
    y <- strsplit(x, "\t")
    names(y) <- sapply(y, `[[`, 1)
    annoList <- lapply(y, `[`, c(-1,-2))
}
