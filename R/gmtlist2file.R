#' write gmt list into gmt file
#' 
#' write gmt list into gmt file
#' 
#' @param gmtlist a list containing gmt
#' @param filename output filename
#' @export
#' 
#' @seealso gmt2list
#' @examples 
#' anno <- "c2.cp.kegg.v5.0.symbols.gmt.xz"
#' annofile <- system.file("extdata", anno, package="cogena")
#' gl <- gmt2list(annofile)
#' gmtfile <- gmtlist2file(gl, filename="")
#' 
gmtlist2file <- function(gmtlist, filename){
    gsname=names(gmtlist) 
    for (i in seq_along(gmtlist) ){ 
        cat(gsname[i], length(gmtlist[[i]]), gmtlist[[i]], file=filename, append=TRUE, sep = "\t")
        cat("\n", append=TRUE, file=filename)
    }
}