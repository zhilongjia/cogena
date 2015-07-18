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
#' @param add2 add 2 clusters information.
#' @details
#' orderMethod:
#' \itemize{
#' \item max. ordered by the max value in clusters beside all
#' \item mean. ordered by the mean value in clusters beside all
#' \item All. ordered by all genes
#' \item I. ordered by the I cluster in two clusters (Up or Down-regulated, add2 should be TRUE)
#' \item II. ordered by the II cluster in two clusters (Up or Down-regulated, add2 should be TRUE)
#' \item a number. like 2, "3".
#' }
#' 
#' @export
#' @docType methods
#' @rdname enrichment
#' @return a matrix with clusters in row and gene-sets in column.
#' @examples
#' data(PD)
#' annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt", 
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
    function(object, method=clusterMethods(object), 
        nCluster=nClusters(object), 
        CutoffNumGeneset=Inf, CutoffPVal=0.05,
        orderMethod="max", roundvalue=TRUE, add2=FALSE) {

    method <- match.arg(method, clusterMethods(object))
    nCluster <- match.arg(nCluster, as.character(nClusters(object)))

    score1 <- object@measures[[method]][[nCluster]]
    score2 <- object@measures[[method]][["2"]]


    if (is.logical(score1)) {
        #warning(paste("For", method, ", the number of clusters:", 
        #nCluster, "Nonexists!"))
        return (score1)
    }

    if (is.logical(score2)) {
        #warning(paste("For", method, ", the number of clusters: 
        #2 Nonexists!"))
        return (score1)
    }

    score1 <- score1[c(as.character(1:nCluster), "All"),]
    score2 <- score2[c("1", "2", "All"),]
    rownames(score2) <- c("I", "II", "All")

    if (nCluster != 2){
        if (isTRUE(add2) ) {
            score <- rbind(score1[1:as.numeric(nCluster),], score2[1:2,], 
                           score1["All",])
            rownames(score) <- c(as.character(1:nCluster), "I", "II", "All")
        } else {
            score <- score1[c(1:as.numeric(nCluster), "All"),]
        }
        
    } else {
        score <- score2
        rownames(score) <- c("I", "II", "All")
    }

    colnames(score) <- tolower(colnames(score))


    # Gene numbers in each clusters, I, II and All.
    NumGeneInCluster <- as.vector( table(geneclusters(object, method, nCluster)) )
    cluster2 <- geneclusters(object, method, "2")
    cluster2_all <- c(length(which(cluster2 == 1)), length(which(cluster2 == 2)), 
            length(geneclusters(object, method, nCluster)) )
    if (nCluster != 2){
        if (isTRUE(add2)) {
            NumGeneInCluster <- c(NumGeneInCluster, cluster2_all)
        }
        
    } else {
        NumGeneInCluster <- cluster2_all
    }
    

    # the orderMethod options
    orderMethod <- match.arg(orderMethod, c(rownames(score), "max", "mean"))
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

    #drop para used as length(index_above_cutoffPVal)==1.
    score <- score[,ncol(score):1, drop=FALSE]

    #order the rownames of score
    #suppressWarnings(score <- 
    #score[order(as.numeric(rownames(score))),,drop=FALSE])
    if (nCluster !=2){
        if (isTRUE(add2)) {
            suppressWarnings(score <- score[c(as.character(1:nCluster), "I", 
                                              "II", "All"),,drop=FALSE])
        } else {
            suppressWarnings(score <- score[c(as.character(1:nCluster), "All"),,drop=FALSE])
        }
        
    } else {
        suppressWarnings(score <- score[c("I", "II", "All"),,drop=FALSE])
    }

    colnames(score) <- tolower(strtrim(colnames(score), 60))
    rownames(score) <- paste(rownames(score), 
        as.character(NumGeneInCluster), sep="#")

    if (roundvalue){
        #-log2(0.5)=1
        #score <- ifelse(round(score,1) < 1, NA, round(score,1))
        score <- round(score,1)
    }
    return (score)
    })
