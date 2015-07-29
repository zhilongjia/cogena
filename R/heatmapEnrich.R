#' heatmap designed for for CMap gene set only
#' 
#' heatmapEnrich is desgined for the cogena result from CMap
#' only so as to collapse the multi-isntance drugs in CMap! 
#' 
#' @inheritParams enrichment
#' @param nCluster as nClust in cogena function.
#' @param orderMethod the order method, default is max, other options are 
#' "mean", "all", "I", "II" or a number meaning the ith cluster.
#' @param CutoffNumGeneset the cut-off of the number of gene sets in the 
#' return table. The default is 20.
#' @param CutoffPVal the cut-off of p-value. The default is 0.05.
#' @inheritParams heatmapPEI
#' @param maintitle a character. Default is "cogena"
#' @param printGS print the enriched gene set names or not. Default is TRUE.
#' 
#' @return a gene set enrichment heatmap
#' 
#' orderMethod:
#' \itemize{
#' \item max. ordered by the max value in clusters beside all
#' \item mean. ordered by the mean value in clusters beside all
#' \item All. ordered by all genes
#' \item Up. ordered by up-regulated genes (add2 should be TRUE)
#' \item Down. ordered by down-regulated genes (add2 should be TRUE)
#' }
#' 
#' @examples
#' \dontrun{
#' heatmapEnrich(clen_res, "kmeans", "3")
#' }
#' 
#' @export
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @docType methods
#' @rdname heatmapEnrich
#' 
setGeneric("heatmapEnrich", 
           function(object, method=clusterMethods(object), 
                    nCluster=nClusters(object), orderMethod="max",
                    CutoffNumGeneset=20, CutoffPVal=0.05, 
                    low="grey", high="red", na.value="white", maintitle="cogena",
                    printGS=TRUE, add2=TRUE) 
               standardGeneric("heatmapEnrich"))

#' @rdname heatmapEnrich
#' @aliases heatmapEnrich,cogena
setMethod("heatmapEnrich", signature(object="cogena"),
          function(object, method=clusterMethods(object), 
                   nCluster=nClusters(object), orderMethod="max",
                   CutoffNumGeneset=20, CutoffPVal=0.05, 
                   low="grey", high="red", na.value="white", maintitle="cogena",
                   printGS=TRUE, add2=TRUE) {
              
              enrichment_score <- enrichment(object, method, nCluster, add2=add2, 
                                             orderMethod="max", CutoffNumGeneset=20, CutoffPVal=0.05)
              
              enrichment_score <- as.data.frame(t(enrichment_score))
              enrichment_score$name <- sapply(strsplit(rownames(enrichment_score) , "_"), "[[", 1)
              meanX <- function(x) {
                  instanceCount <- length(which(x>=round(-log2(CutoffPVal))))
                  meanScore <- mean(x[which(x>=round(-log2(CutoffPVal)))])
                  if (is.nan(meanScore) || instanceCount == 1) {
                      return (0)
                  } else {
                      return (round(meanScore, 1)) 
                  }
              }
              name <- NULL
              score <- dplyr::group_by(enrichment_score, name)
              score <- summarise_each(score, funs(meanX))
              score <- score[which(rowSums(score[,-1])!=0), ]
              rownames(score) <- score$name
              if (nrow(score)==0){
                  stop("No enrichment for this cluster!")
              }
              
              score <- t(subset(as.data.frame(score), select=-name))
              score[score==0] <- NA
              
              # the orderMethod options
              orderMethod <- match.arg(orderMethod, c(rownames(score), "max", "mean"))
              if (orderMethod == "mean") {
                  score = score[,order(colMeans(score, na.rm=TRUE), decreasing=TRUE), drop=FALSE]
                  index_above_cutoffPVal <- which(suppressWarnings(
                      apply(score, 2, max, na.rm=TRUE)) > -log2(CutoffPVal))
              } else if (orderMethod == "max") {
                  colMax <- function(X) {suppressWarnings(
                      apply(X, 2, max, na.rm=TRUE))}
                  score = score[,order(colMax(score), decreasing=TRUE), drop=FALSE]
                  index_above_cutoffPVal <- which(suppressWarnings(
                      apply(score, 2, max, na.rm=TRUE)) > -log2(CutoffPVal))
              } else if (orderMethod %in% rownames(score)) {
                  score = score[, order(score[orderMethod,], decreasing=TRUE), drop=FALSE]
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
              
              score <- score[,ncol(score):1, drop=FALSE]
              if (printGS==TRUE) {
                  cat (rev(colnames(score)), sep ="\t")
              }
              
              cl_color0 <- upDownGene(object, method, nCluster, add2)
              cl_color <- sapply(cl_color0, function(x) {if (x>=0.6) "red" else if (x<=-0.6) "green4" else "black"})
              cl_color <- c(cl_color, "blue")

              enrich_score <- reshape2::melt(score)
              #legend breaks
              if (max(enrich_score$value, na.rm=TRUE) > 15 ){
                  breaks <- seq(15, max(enrich_score$value, na.rm=TRUE), 10)
                  breaks <- c(4.32, breaks)
              } else {
                  breaks <- NULL
              }
              Var1=Var2=value=NULL

              ggplot(enrich_score, aes(as.factor(Var1), Var2)) + 
                  geom_tile(aes(fill = value)) + 
                  scale_fill_gradient2("score",  mid=low, midpoint=4, low=low, 
                                       high=high, na.value=na.value, breaks=breaks) +
                  geom_text(aes(fill=value, label=value),size=4, na.rm=TRUE) +
                  labs(list(title = maintitle, x = "Cluster", y = "Gene set")) +
                  theme(axis.text.y = element_text(size = rel(1.5), face="bold")) +
                  theme(axis.text.x = element_text(size = rel(1.3), angle=30, 
                                                   face="bold", color=cl_color)) 
          }
)

