

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
