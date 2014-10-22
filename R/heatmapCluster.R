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


# @rdname heatmapCluster
#' @exportMethod heatmapCluster
#' @aliases heatmapCluster
setMethod("heatmapCluster", signature(object="clena"),
          function (object, method=clusterMethods(object), nClusters=nClusters(object),
                    sampleColor=c("darkblue", "cyan"), 
                    clusterColor, ...){

              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))
              ColSideColors <- map2col(as.numeric(as.factor(object@sampleLabel)), sampleColor)
              mat <- mat(object)
              colnames(mat) <- paste(colnames(mat), object@sampleLabel, sep="_")
              clusterColor=rainbow(nClusters)
              
              if (method %in% c("hierarchical", "diana", "agnes")) {
                  print (table(cutree(clusters(object, method), nClusters)))
                  mat <- mat[order(cutree(clusters(object, method), nClusters), decreasing=TRUE),]
                  heatmap(mat, Rowv=NA, Colv=NA, labRow=NA,
                          ColSideColors=ColSideColors, ylab="Genes", 
                          main=paste(method, nClusters, "clusters",sep="_"),
                          RowSideColors=map2col(sort(cutree(clusters(object, method), nClusters), decreasing=TRUE), clusterColor),
                          ...)
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
                          RowSideColors=map2col(sort(clusters(object, method)[[nClusters]]$cluster, decreasing=TRUE), clusterColor),
                          ...)
              }
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