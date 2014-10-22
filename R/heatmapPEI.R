#' heatmap.PEI
#'
#' heatmap.PEI
#' @param object a clean object
#' @param orderMethod the criterion of ordering the heatmap. The default is max, the other options is mean
#' @method heatmap.PEI
#' @export heatmap.PEI
setMethod("heatmap.PEI", signature(object="clena"),
          function(object, method=clusterMethods(object), 
                   nClusters=nClusters(object), CutoffNumGeneset=25, orderMethod="max") {
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))
              enrichment <- enrichment(object, method, nClusters)
              
              #annotate the cluster with number of genes
              #NumGeneInCluster <- vector(mode="character",length=nClusters+1)
              if (method %in% c("hierarchical", "diana", "agnes")) {
                  NumGeneInCluster <- as.vector(table(cutree(clusters(object, method), k=nClusters)))
                  NumGeneInCluster <- c(NumGeneInCluster, length(cutree(clusters(object, method), k=nClusters)))
              } else 
              { NumGeneInCluster <- as.vector(table(clusters(object, method)[[nClusters]]$cluster))
                NumGeneInCluster <- c(NumGeneInCluster, length(clusters(object, method)[[nClusters]]$cluster))
              }
              cat ("#Gene In each Cluster:")
              cat (NumGeneInCluster)
              
              # the orderMethod options to order the enrichment: mean and max
              if (orderMethod == "mean") {
                  enrichment = enrichment[,order(colMeans(enrichment[-nrow(enrichment),], na.rm=TRUE), decreasing=TRUE)]
              } else if (orderMethod == "max") {
                  colMax <- function(X) {suppressWarnings(apply(X, 2, max, na.rm=TRUE))}
                  enrichment = enrichment[,order(colMax(enrichment[-nrow(enrichment),]), decreasing=TRUE)]
              } else {
                  warning(paste("\n wrong orderMethod:", orderMethod)); break
              }
              #Trim the NO. of Geneset no more than CutoffNumGeneset
              if (ncol(enrichment)>CutoffNumGeneset) {enrichment = enrichment[,1:CutoffNumGeneset]}
              
              # if (ncol(enrichment)>CutoffNumGeneset) {
              # 	colMax <- function(X) {apply(X, 2, max, na.rm=TRUE)}
              # 	enrichment = enrichment[,order(colMax(enrichment[-nrow(enrichment),]), decreasing=TRUE)[1:CutoffNumGeneset]]
              # 	#enrichment = enrichment[,order(colMeans(enrichment[-nrow(enrichment),], na.rm=TRUE), decreasing=TRUE)[1:CutoffNumGeneset]]
              # }
              enrichment <- enrichment[,ncol(enrichment):1]
              colnames(enrichment) <- tolower(strtrim(colnames(enrichment), 60))
              rownames(enrichment) <- paste(rownames(enrichment), as.character(NumGeneInCluster), sep="#")
              enrichment <- ifelse(round(enrichment,1)==0, NA, round(enrichment,1))
              enrichment <- reshape2::melt(enrichment)
              ggplot(enrichment, aes(as.factor(Var1), Var2)) + 
                  geom_tile(aes(fill = value)) + 
                  scale_fill_gradient(low="green", high="red", na.value="white") +
                  geom_text(aes(fill=value, label=value),size=4, na.rm=TRUE) +
                  labs(list(title = paste(method, nClusters,"Clusters-Set-Analysis"), x = "Cluster", y = "Set")) +
                  theme(axis.text.y = element_text(size = rel(1.4), face="bold")) +
                  theme(axis.text.x = element_text(size = rel(1.3), angle=30)) 
          })
