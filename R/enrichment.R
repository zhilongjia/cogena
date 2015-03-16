#' get the enrichment table from a cogena object.
#' 
#' get the enrichment table from a cogena object with certain clustering 
#' methods and number of clusters.
#' 
#' @inheritParams clusters
#' @param nClusters as nClust in cogena function.
#' @param CutoffNumGeneset the cut-off of the number of gene sets in the 
#' return table
#' @param CutoffPVal the cut-off of p-value. The default is 0.05.
#' @param orderMethod the order method, default is max, other options are 
#' "mean", "all", "I", "II"
#' @param roundvalue The default is TRUE. whether or not round the data. 
#' such as round(1.54, 1)=1.5
#' 
#' @details
#' orderMethod:
#' \itemize{
#' \item max. ordered by the max value in clusters beside all
#' \item mean. ordered by the mean value in clusters beside all
#' \item all. ordered by all genes
#' \item I. ordered by the I cluster in only two clusters (Up or Down-regulated)
#' \item II. ordered by the II cluster in only two clusters (Up or Down-regulated)
#' }
#' 
#' @export
#' @docType methods
#' @rdname enrichment
#' @return a matrix with clusters in row and gene-sets in column.
#' @examples
#' data(PD)
#' enrichment.table1 <- enrichment(cogena_result, "kmeans", "3")
#' enrichment.table2 <- enrichment(cogena_result, "kmeans", "3", 
#' CutoffNumGeneset=10, orderMethod="mean")
setGeneric("enrichment", function(object, method, nClusters, 
                                  CutoffNumGeneset=Inf, 
                                  CutoffPVal=0.05, orderMethod="max", 
                                  roundvalue=TRUE) 
    standardGeneric("enrichment"))

#' @rdname enrichment
#' @aliases enrichment,cogena_methods
setMethod("enrichment", signature(object="cogena"),
          function(object, method=clusterMethods(object), 
                   nClusters=nClusters(object), 
                   CutoffNumGeneset=Inf, CutoffPVal=0.05,
                   orderMethod="max", roundvalue=TRUE) {
              
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, 
                                     as.character(nClusters(object)))
              score1 <- object@measures[[method]][[nClusters]]
              score2 <- object@measures[[method]][["2"]]

              if (is.logical(score1)){
                  #warning(paste("For", method, ", the number of clusters:", nClusters, "Nonexists!"))
                  return (score1)
              }

              if (is.logical(score2)){
                  #warning(paste("For", method, ", the number of clusters: 
                  #2 Nonexists!"))
                  return (score1)
              }

              if (nClusters != 2){
                  score <- rbind(score1[1:as.numeric(nClusters),], score2[1:2,]
                                 , score1["All",])
                  rownames(score) <- c(as.character(1:nClusters), "I", "II", 
                                       "All")
              } else {
                  score <- score2
                  rownames(score) <- c("I", "II", "All")
              }

              colnames(score) <- tolower(colnames(score))


              if (method %in% c("hierarchical", "diana", "agnes")) {
                  NumGeneInCluster <- as.vector(
                      table(cutree(clusters(object, method), k=nClusters)))
                  cluster2 <- cutree(clusters(object, method), 2)
                  cluster2_all <- c(length(which(cluster2 == 1)), 
                                    length(which(cluster2 == 2)), 
                                    length(cutree(clusters(object, method), 
                                                  k=nClusters)))
                  if (nClusters != 2) {
                      NumGeneInCluster <- c(NumGeneInCluster, cluster2_all)
                  } else {
                      NumGeneInCluster <- cluster2_all
                  }
              } else 
              { NumGeneInCluster <- as.vector(
                  table(clusters(object, method)[[nClusters]]$cluster))
                cluster2 <- clusters(object, method)[["2"]]$cluster
                cluster2_all <- c(length(which(cluster2 == 1)), 
                                  length(which(cluster2 == 2)), 
                                  length(clusters(object, method)[[nClusters]]$cluster))
                if (nClusters != 2){
                    NumGeneInCluster <- c(NumGeneInCluster, cluster2_all)
                } else {
                    NumGeneInCluster <- cluster2_all
                }
              }

              # the orderMethod options to order the score: mean, all and max
              if (orderMethod == "mean") {
                  score = score[,order(colMeans(score, na.rm=TRUE), decreasing=TRUE)]
              } else if (orderMethod == "max") {
                  colMax <- function(X) {suppressWarnings(apply(X, 2, max, na.rm=TRUE))}
                  score = score[,order(colMax(score), decreasing=TRUE)]
                  #score = score[,order(colMax(score[-nrow(score),]), decreasing=TRUE)]
              } else if (orderMethod == "all") {
                  score = score[, order(score["All",], decreasing=TRUE)]
              } else if (orderMethod == "I") {
                  score = score[, order(score["I",], decreasing=TRUE)]
              } else if (orderMethod == "II") {
                  score = score[, order(score["II",], decreasing=TRUE)]
              } else {
                  warning(paste("\n wrong orderMethod:", orderMethod)); break
              }

              #Trim the NO. of Geneset no more than CutoffNumGeneset and delete the all NAs
              index_above_cutoffPVal <- which(suppressWarnings(apply(score, 2, max, na.rm=TRUE)) > -log2(CutoffPVal))
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
              #suppressWarnings(score <- score[order(as.numeric(rownames(score))),,drop=FALSE])
              if (nClusters !=2){
                  suppressWarnings(score <- score[c(as.character(1:nClusters), "I", "II", "All"),,drop=FALSE])
              } else {
                  suppressWarnings(score <- score[c("I", "II", "All"),,drop=FALSE])
              }

              colnames(score) <- tolower(strtrim(colnames(score), 60))
              rownames(score) <- paste(rownames(score), as.character(NumGeneInCluster), sep="#")

              if (roundvalue){
                  #-log2(0.5)=1
                  #score <- ifelse(round(score,1) < 1, NA, round(score,1))
                  score <- round(score,1)
              }

              return (score)
          })
