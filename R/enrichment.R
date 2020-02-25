#' get the enrichment table from a cogena object.
#' 
#' get the enrichment table from a cogena object with certain clustering 
#' methods and number of clusters.
#' 
#' @inheritParams geneclusters
#' @param nCluster as nClust in cogena function.
#' @param CutoffNumGeneset the cut-off of the number of gene sets in the 
#' return table
#' @param CutoffPVal the cut-off of p-value. The default is 0.05.
#' @param orderMethod the order method, default is max, other options are 
#' "mean", "all", "I", "II" or a number meaning the ith cluster.
#' @param roundvalue The default is TRUE. whether or not round the data. 
#' such as round(1.54, 1)=1.5
#' @param add2 enrichment score for add Up and Down reuglated genes.
#' @details
#' orderMethod:
#' \itemize{
#' \item max. ordered by the max value in clusters beside all
#' \item mean. ordered by the mean value in clusters beside all
#' \item All. ordered by all genes
#' \item I. ordered by the I cluster in two clusters (Up or Down-regulated, add2 should be TRUE)
#' \item II. ordered by the II cluster in two clusters (Up or Down-regulated, add2 should be TRUE)
#' \item a character number. like "3".
#' }
#' 
#' @export
#' @docType methods
#' @rdname enrichment
#' @return a matrix with clusters in row and gene-sets in column.
#' @examples
#' data(Psoriasis)
#' annofile <- system.file("extdata", "c2.cp.kegg.v7.01.symbols.gmt.xz", 
#' package="cogena")
#' 
#' \dontrun{
#' genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
#' 
#' enrichment.table1 <- enrichment(clen_res, "kmeans", "3")
#' enrichment.table2 <- enrichment(clen_res, "kmeans", "3", 
#' CutoffNumGeneset=10, orderMethod="mean")
#' }
#' 
setGeneric("enrichment", function(object, method, nCluster, 
    CutoffNumGeneset=Inf, CutoffPVal=0.05, orderMethod="max", roundvalue=TRUE,
    add2=FALSE) 
    standardGeneric("enrichment"))

#' @rdname enrichment
#' @aliases enrichment,cogena_methods
setMethod("enrichment", signature(object="cogena"),
    function(object, method, nCluster, 
        CutoffNumGeneset=Inf, CutoffPVal=0.05,
        orderMethod="max", roundvalue=TRUE, add2=TRUE) {

    method <- match.arg(method, clusterMethods(object))
    nCluster <- match.arg(nCluster, as.character(nClusters(object)))

    score1 <- object@measures[[method]][[nCluster]]
    
    if (is.logical(score1) ) {
        # warning(paste("For", method, ", the number of clusters:", nCluster, "Nonexists!"))
        return (score1)
    } else if (is.null(score1)) {
        stop(paste("For", method, "with ", nCluster, "clusters: Nonexists in the input cogena object!"))
    }

    # Gene numbers in each clusters, I, II and All.
    NumGeneInCluster <- as.vector( table(geneclusters(object, method, nCluster)) )
    NumGeneIncluster2 <- c(length(object@upDn$upGene), length(object@upDn$dnGene))
    
    if (isTRUE(add2)) {
        score <- score1[c(1:as.numeric(nCluster), c("Up", "Down", "All")),]
        NumGeneInCluster <- c(NumGeneInCluster, NumGeneIncluster2, length(geneclusters(object, method, nCluster)))
    } else {
        score <- score1[c(1:as.numeric(nCluster), "All"),]
        NumGeneInCluster <- c(NumGeneInCluster, length(geneclusters(object, method, nCluster)))
    }
    # colnames(score) <- tolower(colnames(score))
    

    # the orderMethod options
    orderMethod <- match.arg(orderMethod, c(rownames(score), "max", "mean", "each"))
    if (orderMethod == "mean") {
        score = score[,order(colMeans(score, na.rm=TRUE), decreasing=TRUE)]
        index_above_cutoffPVal <- which(suppressWarnings(
            apply(score, 2, max, na.rm=TRUE)) > -log2(CutoffPVal))
    } else if (orderMethod == "max") {
        colMax <- function(X) {suppressWarnings(
        apply(X, 2, max, na.rm=TRUE))}
            score = score[,order(colMax(score), decreasing=TRUE)]
        index_above_cutoffPVal <- which(suppressWarnings(
            apply(score, 2, max, na.rm=TRUE)) > -log2(CutoffPVal))
    } else if (orderMethod %in% rownames(score)) {
        score = score[, order(score[orderMethod,], decreasing=TRUE)]
        index_above_cutoffPVal <- 
        which(score[orderMethod,] > -log2(CutoffPVal))
    } else if (orderMethod == "each") {
        gs_per_cluster <- round(CutoffNumGeneset/as.numeric(nCluster))
        gs_i <- list()
        for (i in rownames(score)) {
            # print (i)
            gs_i[[i]] <- names(sort(score[i,], decreasing = TRUE))[1:gs_per_cluster]
        }
        score <- score[,unlist(gs_i)]
        index_above_cutoffPVal <- which(suppressWarnings(
            apply(score, 2, max, na.rm=TRUE)) > -log2(CutoffPVal))
    }

    if (length(index_above_cutoffPVal) > CutoffNumGeneset){
        score <- score[,c(1:CutoffNumGeneset)]
    } else if (length(index_above_cutoffPVal) == 0){
        score <- NA
        return (score)
    } else {
        #drop para used as length(index_above_cutoffPVal)==1.
        score <- score[,index_above_cutoffPVal, drop=FALSE]
    }

    # drop para used as length(index_above_cutoffPVal)==1.
    score <- score[,ncol(score):1, drop=FALSE]

    # Upper cell type and conc in CMAP
    # if (grepl("@", colnames(score)[1])) {
    #     colnames(score) <- colnames(score)
    #     # colnames(score) <- paste(sapply(strsplit(colnames(score), "@"), "[", 1), 
    #     #                          toupper(sapply(strsplit(colnames(score), "@"), "[", 2)), 
    #     #                          sep="@")
    # } else {
    #     # colnames(score) <- tolower(strtrim(colnames(score), 60))
    #     colnames(score) <- tolower(colnames(score))
    # }
    
    rownames(score) <- paste(rownames(score), as.character(NumGeneInCluster), sep="#")
    
    if (roundvalue){
        score <- round(score,1)
    }
    
    return (score)

})
