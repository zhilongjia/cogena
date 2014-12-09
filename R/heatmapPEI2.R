#' heatmap of the gene set enrichment_score matrix directly
#' 
#' heatmap of the gene set enrichment_score matrix directly. After obtaining the 
#' ennrichemt of clusters in the gene sets via \code{\link{enrichment}}, the
#' heatmapPEI2 will show it as a heatmap.
#' 
#' @inheritParams heatmapPEI
#' @whichCluster which cluster should be based to filter, value like "2#68".
#' 
#' @details
#' This function aims to heatmap the enrichment_score directly. This is helpful
#' on condition that there are so many enriched gene sets and you can filter the
#' enrichment_score based on a criteria, like just one cluster. 
#' @examples
#' ##
#' @export
#' @docType methods
#' @rdname heatmapPEI2
setGeneric("heatmapPEI2", function(object, enrichment_score, ...) standardGeneric("heatmapPEI2"))


#' @aliases heatmapPEI2,clena
setMethod("heatmapPEI2", signature(object="clena"),
          function(object, enrichment_score, method, nClusters, whichCluster,
                   CutoffNumGeneset=60, 
                   #CutoffPVal=0.05,
                   #orderMethod="max", roundvalue=TRUE,
                   low="grey", high="red", na.value="white") {
              
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))
              whichCluster <- match.arg(whichCluster, as.character(rownames(enrichment_score)))

              enrichment_score <- enrichment_score[,which(!is.na(enrichment_score[whichCluster,]))]
              enrichment_score <- enrichment_score[, order(enrichment_score[whichCluster,], decreasing=FALSE)]
              if (ncol(enrichment_score)> CutoffNumGeneset){
                  enrichment_score <- enrichment_score[, (ncol(enrichment_score)-CutoffNumGeneset):ncol(enrichment_score)]
              }

              if (ncol(enrichment_score)==0){
              	stop("No enrichment for this cluster!")
              }

              enrichment_score <- reshape2::melt(enrichment_score)
              #legend breaks
              if (max(enrichment_score$value, na.rm=TRUE) > 15 ){
                  breaks <- seq(15, max(enrichment_score$value, na.rm=TRUE), 10)
                  breaks <- c(4.32, breaks)
              } else {
                  breaks <- NULL
                  
              }
              ggplot(enrichment_score, aes(as.factor(Var1), Var2)) + 
                  geom_tile(aes(fill = value)) + 
                  scale_fill_gradient2("score",  mid=low, midpoint=4, low=low, high=high, na.value=na.value, breaks=breaks) +
                  geom_text(aes(fill=value, label=value),size=4, na.rm=TRUE) +
                  labs(list(title = paste(method, nClusters,"Clusters-Set-Analysis"), x = "Cluster", y = "Set")) +
                  theme(axis.text.y = element_text(size = rel(1.5), face="bold")) +
                  theme(axis.text.x = element_text(size = rel(1.3), angle=30, face="bold")) 
          })
