#' read gmt file as a list
#' 
#' read Gene Matrix Transposed (gmt) file and output a list with the the first 
#' column as the names of items in the list. see 
#' \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{Gene Matrix Transposed file format}
#' for more details.
#' 
#' @param annofile a gmt file. Examples are from MSigDB Collections.
#' A list of gene set could be find in the vignette of cogena
#' @return a gmt list
#' @export
#' 
#' @seealso gmtlist2file
#' @examples
#' anno <- "c2.cp.kegg.v5.0.symbols.gmt.xz"
#' annofile <- system.file("extdata", anno, package="cogena")
#' gl <- gmt2list(annofile)
#' 
gmt2list <- function(annofile){
    if (!file.exists(annofile)) {
        stop("There is no such gmt file.")
    }
    
    if (tools::file_ext(annofile) == "xz") {
        annofile <- xzfile(annofile)
        x <- scan(annofile, what="", sep="\n", quiet=TRUE)
        close(annofile)
    } else if (tools::file_ext(annofile) == "gmt") {
        x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    } else {
        stop ("Only gmt and gmt.xz are accepted for gmt2list")
    }
    
    y <- strsplit(x, "\t")
    names(y) <- sapply(y, `[[`, 1)
    
    annoList <- lapply(y, `[`, c(-1,-2))
}

