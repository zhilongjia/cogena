#' Basic methods for a cogena object.
#' 
#' clusterMethods: get the methods of clustering used.
#' @param object a cogena object
#' @exportMethod clusterMethods
#' @import methods
#' @docType methods
#' @rdname cogena_methods
#' @return clusterMethods: a character vector.
#' @examples
#' data(PD)
#' annofile <- system.file("extdata", "c2.cp.kegg.v4.0.symbols.gmt", 
#' package="cogena")
#' cogena_result <- cogena(DEexprs, nClust=2:3, 
#' clMethods=c("hierarchical","kmeans"), metric="correlation", 
#' method="complete",  annofile=annofile, sampleLabel=sampleLabel, 
#' ncore=1, verbose=TRUE)
#' clusterMethods(cogena_result)
#' 
#' 
setGeneric("clusterMethods", function(object) standardGeneric("clusterMethods"))

#' @rdname cogena_methods
#' @aliases clusterMethods,cogena_methods
setMethod("clusterMethods",signature(object="cogena"),
          function(object) return(object@clMethods))


#' nClusters: get the number of clusters from a cogena object.
#' @inheritParams clusterMethods
#' @examples
#' \dontrun{
#' nClusters(cogena_result)
#' }
#' @return nClusters: a numeric vector.
#' @export nClusters
#' @docType methods
#' @rdname cogena_methods
setGeneric("nClusters", function(object) standardGeneric("nClusters"))

#' @rdname cogena_methods
#' @aliases nClusters,cogena_methods
setMethod("nClusters",signature(object="cogena"),
          function(object) return(object@nClust))



#' clusters: get the cluster of a certain clustering method.
#' @inheritParams clusterMethods
#' @param method as clMethods in cogena function
#' @return clusters: a list or hclust depends on the method
#' @examples
#' \dontrun{
#' clusters(cogena_result, "kmeans")
#' clusters(cogena_result, "hierarchical")
#' }
#' @export clusters
#' @docType methods
#' @rdname cogena_methods
setGeneric("clusters", function(object, method) standardGeneric("clusters"))

#' @rdname cogena_methods
#' @aliases clusters,cogena_methods
setMethod("clusters",signature(object="cogena"),
          function(object, method=clusterMethods(object)) {
              method <- match.arg(method, clusterMethods(object))
              return(object@clusterObjs[[method]])})


#' mat: get the original data from a cogena object.
#' @inheritParams clusterMethods
#' @examples
#' \dontrun{
#' mat(cogena_result)
#' }
#' @return mat: a matrix
#' @export mat
#' @docType methods
#' @rdname cogena_methods
setGeneric("mat", function(object) standardGeneric("mat"))

#' @rdname cogena_methods
#' @aliases mat,cogena_methods
setMethod("mat",signature(object="cogena"),
          function(object) return(object@mat))


#' summary: a summary of a cogena object.
#' @return summary: a summary of a cogena object.
#' @examples
#' \dontrun{
#' summary(cogena_result)
#' }
#' @rdname cogena_methods
#' @exportMethod summary
setMethod("summary","cogena",
          function(object) {
              cat("\nClustering Methods:\n",clusterMethods(object),"\n\n")
              cat("The Number of Clusters:\n",nClusters(object),"\n\n")
              cat("Metric of Distance Matrix:\n", object@metric, "\n\n")
              cat("Agglomeration method for hierarchical clustering (hclust and agnes):\n", object@method, "\n\n")
          })
