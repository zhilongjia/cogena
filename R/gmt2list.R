#' read gmt file as a list
#' 
#' read Gene MatriX (gmt) file and output a list with the the first column as 
#' the names of items in the list.
#' 
#' @param anno a gmt file. Examples are from MSigDB Collections.
#' The default is c2.cp.v4.0.symbols.gmt. Others like c2.cgp.v4.0.symbols.gmt, 
#' c2.cp.kegg.v4.0.symbols.gmt, c2.cp.reactome.v4.0.symbols.gmt, 
#' c5.bp.v4.0.symbols.gmt, c6.all.v4.0.symbols.gmt. Or a gmt format file.
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
