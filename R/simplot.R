#' plotting similarity matrix
#' 
#' @inheritParams geneclusters
#' @param nCluster as nClust in cogena function.
#' @param CutoffNumGeneset the cut-off of the number of gene sets in the 
#' return table
#' @param decreasing parameter for order
#' @import corrplot
#' @export simplot
#' @examples 
#' 
#' 
#' genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#'     
#' annoGMT <- "c2.cp.kegg.v5.0.symbols.gmt.xz"
#' annofile <- system.file("extdata", annoGMT, package="cogena")
#' 
#' data(Psoriasis)
#' clMethods <- c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes")

#' genecl_result <- coExp(DEexprs, nClust=10, clMethods="pam", 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)   

#' cogena_res <- clPage(genecl_result, "pam", "10", sampleLabel=sampleLabel)
#' simplot(cogena_res, "pam", "10")
simplot <- function(object, method, nCluster, orderMethod="min", CutoffNumGeneset=20, decreasing=FALSE){
    
    sim <- object@annotation
    score <- t(sim)
    rownames(score) <- sub("#.*", "", colnames(sim))
    
    orderMethod <- match.arg(orderMethod, c(rownames(score), "max", "mean", "min"))
    
    if (orderMethod == "mean") {
        score1 = score[,order(colMeans(score, na.rm=TRUE), decreasing=decreasing)]
    } else if (orderMethod == "max") {
        colMax <- function(X) {suppressWarnings(
            apply(X, 2, max, na.rm=TRUE))}
        score1 = score[,order(colMax(score), decreasing=TRUE)]
    } else if (orderMethod %in% rownames(score)) {
        score1 = score[, order(score[orderMethod,], decreasing=decreasing)]
    } else if (orderMethod == "min") {
        colMin <- function(X) {suppressWarnings(
            apply(X, 2, min, na.rm=TRUE))}
        score1 = score[,order(colMin(score), decreasing=FALSE)]
    }
    
    score2 <- score1[, 1:CutoffNumGeneset]
    rownames(score2) <- colnames(sim)
    
    cl_color0 <- upDownGene(object, method, nCluster, add2=TRUE)
    cl_color <- sapply(cl_color0, function(x) {if (x>=0.6) "red" else if (x<=-0.6) "green4" else "black"})
    cl_color <- c(cl_color, "blue")
    
    corrplot::corrplot(t(score2), is.corr=FALSE,  tl.x.col=cl_color, cl.lim=c(-1, 1) )
}

