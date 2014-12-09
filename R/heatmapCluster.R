#' heatmap of gene expression profilings with cluster indication.
#' 
#' heatmap of gene expression profilings with cluster-based color indication.
#' @inheritParams enrichment
#' @param sampleColor a color vector with the sample length. The default is c("darkblue", "cyan").
#' @param clusterColor a color vector with the cluster length. The default is rainbow(nClusters(object)).
#' @param ... other parameters to stats::heatmap, such as col=topo.colors(100).
#' @export
#' @docType methods
#' @usage heatmapCluster(object, method, nClusters, sampleColor = c("darkblue", "cyan"), 
#' clusterColor = rainbow(nClusters), ...)
#' @seealso \code{\link{clena}} and \code{\link{heatmapPEI}}
#' @examples 
#' #summay this clena object
#' summary(GSE7621.SAM.clena.cluster)
#'
#' #heatmapCluster
#' heatmapCluster(GSE7621.SAM.clena.cluster, "hierarchical", "12", col=topo.colors(100))
#' 
setGeneric("heatmapCluster", 
           function(object, method, nClusters, sampleColor=c("darkblue", "cyan"), 
                    clusterColor, ...) standardGeneric("heatmapCluster"))

#load("/home/zjia/data/Jobbing/PD/result/PD_clena_analysis.c2.cp.kegg.v4.0.symbols.gmt.RData")

# @rdname heatmapCluster
#' @exportMethod heatmapCluster
#' @aliases heatmapCluster
setMethod("heatmapCluster", signature(object="clena"),
          function (object, method=clusterMethods(object), nClusters=nClusters(object),
                    sampleColor=c("darkblue", "cyan"), 
                    clusterColor, ...){
              
              #get the parameters
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))
              mat <- mat(object)
              
              if (method %in% c("hierarchical", "diana", "agnes")) {
                  cluster_size <- cutree(clusters(object, method), nClusters)
                  names(cluster_size) <- rownames(mat)
                  cluster_size2 <- cutree(clusters(object, method), "2")
                  names(cluster_size2) <- rownames(mat)
              } else {
                  cluster_size <- clusters(object, method)[[nClusters]]$cluster
                  cluster_size2 <- clusters(object, method)[["2"]]$cluster
              }
              print (table(cluster_size2))
              if (nClusters != "2"){
                  print (table(cluster_size))
                }
              #print (cluster_size)
              #reorder the mat based on the clustering and type of sample.
              sampleLabel <- sort(object@sampleLabel)
              mat <- mat[order(cluster_size, decreasing=FALSE), names(sampleLabel)]
              #add the type of sample into the colnames
              #colnames(mat) <- paste(colnames(mat), sampleLabel, sep="_")
              
              #color setting
              ColSideColors <- map2col(as.numeric(as.factor(sampleLabel)), sampleColor)

              clusterColor <- sample(rainbow(nClusters, alpha = c(1, 0.6)))
              clusterColor2 <- c("coral3", "deepskyblue1")#sample(rainbow(2, alpha = c(1, 0.6)))
              #clusterColor <- rainbow(nClusters, alpha = c(1, 0.6))
              #make the neighboor colors different
              #clusterColor <- clusterColor[clusterColor]
              
              RowSideColors <- map2col(sort(cluster_size, decreasing=FALSE), clusterColor)
              RowSideColors2 <- map2col(cluster_size2[names(sort(cluster_size, decreasing=FALSE))], clusterColor2)
              
              if (nClusters != "2"){
                  RowSideColors <- t(cbind(RowSideColors2, RowSideColors))
                  rownames(RowSideColors) <- paste("Size:", c(2, nClusters))
                  } else {
                  	RowSideColors <- t(as.matrix(RowSideColors))
                  }

              heatmap.3(mat, col=redgreen(75) , trace="none", scale="row", 
                        Rowv=FALSE, Colv=FALSE, dendrogram="none", labRow=NA, 
                        colsep=length(which(sampleLabel==sampleLabel[1])), 
                        rowsep=cumsum( table(cluster_size)), #adjCol=c(0.8,0),
                        sepcolor="white", sepwidth=c(0.05,1),
                        key=TRUE, symkey=FALSE, density.info="none", keysize=1.5, 
                        main=paste(method, nClusters, "clusters",sep="_"),
                        ColSideColors=as.matrix(ColSideColors), ylab="Genes", cexCol=0.8,
                        RowSideColors=RowSideColors, RowSideColorsSize=ifelse(nClusters!=2,2,1))
              par(lend = 0, xpd=TRUE)
              legend("left", legend = paste0(1:nClusters),
                         col = clusterColor, lty= 1, lwd = 20, bty = "n", title = paste(nClusters, "Clusters"))
              if (nClusters != "2"){
                  legend("bottomleft", legend = as.character(as.roman(1:2)),
                         col = clusterColor2, lty= 1, lwd = 20, bty = "n", title = paste(2, "Clusters"))
              }
              legend("top", legend = names(table(sampleLabel)), col = sampleColor, 
                     lty=1, lwd=20, bty = "n", title = "Type of Sample")
          })


################################################################################
#map2color: get the color vector from the numeric vector x using pal, such as, rainbow(200)
#Example: map2color(cutree(clusters(GSE48350.clena.cluster, "hierarchical"), 3),rainbow(200))
#http://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r
##############################################################################
map2col<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}