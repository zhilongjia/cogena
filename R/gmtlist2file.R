#' write gmt list into gmt file
#' 
#' write gmt list into gmt file
#' 
#' @param gmtlist a list containing gmt
#' @param filename output filename
#' @param description a detailed description vector of a gene set this line. 
#' length of the vector should be the same as the length of gmtlist.
#' @return NA
#' 
#' @export
#' 
#' @seealso gmt2list
#' @examples 
#' anno <- "c2.cp.kegg.v7.01.symbols.gmt.xz"
#' annofile <- system.file("extdata", anno, package="cogena")
#' gl <- gmt2list(annofile)
#' gmtfile <- gmtlist2file(gl, filename="" )
#' 
gmtlist2file <- function(gmtlist, filename, description=NULL){
    gsname=names(gmtlist) 
    
    if (file.exists(filename) ) { file.remove(filename) }
    
    if (is.null(description)) {
        for (i in seq_along(gmtlist) ){ 
            cat(gsname[i], length(gmtlist[[i]]), gmtlist[[i]], file=filename, append=TRUE, sep = "\t")
            cat("\n", append=TRUE, file=filename)
        }
    } else if (length(description) == length(gmtlist) ) {
        description <- gsub("\r\n|\r|\n","",description)
        for (i in seq_along(gmtlist) ){ 
            cat(gsname[i], description[i], gmtlist[[i]], file=filename, append=TRUE, sep = "\t")
            cat("\n", append=TRUE, file=filename)
        }
    } else {
        stop("The length of description should be the same as the length of gmtlist.")
    }
    
}
