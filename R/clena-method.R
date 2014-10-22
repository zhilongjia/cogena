


#' heatmap.PEI
#'
#' heatmap.PEI
#' @export heatmap.PEI
setGeneric("heatmap.PEI", function(object, ...) standardGeneric("heatmap.PEI"))

#' heatmap.PEI
#'
#' heatmap.PEI
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

##############################################################################
setGeneric("heatmap.Cluster", function(object, method, nClusters, ...) standardGeneric("heatmap.Cluster"))
#' @export heatmap.Cluster
setMethod("heatmap.Cluster", signature(object="clena"),
          function (object, method=clusterMethods(object), nClusters=nClusters(object)){
              ################################################################################
              #map2color: get the color vector from the numeric vector x using pal, such as, rainbow(200)
              #Example: map2color(cutree(clusters(GSE48350.clena.cluster, "hierarchical"), 3),rainbow(200))
              #http://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r
              ##############################################################################
              map2col<-function(x,pal,limits=NULL){
                  if(is.null(limits)) limits=range(x)
                  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
              }
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))
              ColSideColors <- map2col(as.numeric(as.factor(object@sampleLabel)), c("darkblue", "cyan"))
              mat <- mat(object)
              colnames(mat) <- paste(colnames(mat), object@sampleLabel, sep="_")
              
              if (method %in% c("hierarchical", "diana", "agnes")) {
                  print (table(cutree(clusters(object, method), nClusters)))
                  mat <- mat[order(cutree(clusters(object, method), nClusters), decreasing=TRUE),]
                  heatmap(mat, Rowv=NA, Colv=NA, labRow=NA,
                          ColSideColors=ColSideColors, ylab="Genes", 
                          main=paste(method, nClusters, "clusters",sep="_"),
                          RowSideColors=map2col(sort(cutree(clusters(object, method), nClusters), decreasing=TRUE), rainbow(nClusters)))
              } else {
                  print (table(clusters(object, method)[[nClusters]]$cluster))
                  #kck <- clusters(object, method)[[nClusters]]$cluster
                  #print(method); print (kck[head(rownames(mat))]); print (kck[tail(rownames(mat))])
                  #head1 <- head(rownames(mat)); tail1 <- tail(rownames(mat))
                  #mat <- mat[names(sort(clusters(object, method)[[nClusters]]$cluster, decreasing=TRUE)),]
                  mat <- mat[order(clusters(object, method)[[nClusters]]$cluster, decreasing=TRUE),]
                  #mat_old <- mat
                  
                  #kck1 <- order(clusters(object, method)[[nClusters]]$cluster)
                  #mat <- mat[kck1,]
                  #print (kck[head1]); print (kck[tail1])
                  
                  heatmap(mat, Rowv=NA, Colv=NA, labRow=NA,
                          ColSideColors=ColSideColors, ylab="Genes", 
                          main=paste(method, nClusters, "clusters",sep="_"),
                          RowSideColors=map2col(sort(clusters(object, method)[[nClusters]]$cluster, decreasing=TRUE),rainbow(nClusters))
                  )
              }
          })

############################################################################
setGeneric("geneInCluster", function(object, ...) standardGeneric("geneInCluster"))
#' @export geneInCluster
setMethod("geneInCluster", signature(object="clena"),
          function (object, method=clusterMethods(object), nClusters=nClusters(object), ith){
              #ith is the ith cluster enquerying
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))
              ith <- match.arg(ith, as.character(seq(1:ith)))
              ith <- as.numeric(ith)
              if (method %in% c("hierarchical", "diana", "agnes")) {
                  gene <- names(which (cutree(clusters(object, method), k=nClusters)==ith))
              } else {
                  gene <- names(which(clusters(object, method)[[nClusters]]$cluster==ith))
              }
              return (gene)
          })