#' Basic methods for a clena object.
#' 
#' clusterMethods: get the methods of clustering used.
#' @param object a clena object
#' @export clusterMethods
#' @docType methods
#' @rdname clena_methods
#' @return clusterMethods: a character vector.
#' @examples
#' data(PD)
#' clusterMethods(clena_result)
#' 
setGeneric("clusterMethods", function(object) standardGeneric("clusterMethods"))


#' @aliases clusterMethods,clena_methods
setMethod("clusterMethods",signature(object="clena"),
          function(object) return(object@clMethods))


#' nClusters: get the number of clusters from a clena object.
#' @inheritParams clusterMethods
#' @examples
#' nClusters(clena_result)
#' @return nClusters: a numeric vector.
#' @export nClusters
#' @docType methods
#' @rdname clena_methods
setGeneric("nClusters", function(object) standardGeneric("nClusters"))


#' @aliases nClusters,clena_methods
setMethod("nClusters",signature(object="clena"),
          function(object) return(object@nClust))



#' clusters: get the cluster of a certain clustering method.
#' @inheritParams clusterMethods
#' @param method as clMethods in clena function
#' @return clusters: a list or hclust depends on the method
#' @examples
#' clusters(clena_result, "kmeans")
#' clusters(clena_result, "hierarchical")
#' @export clusters
#' @docType methods
#' @rdname clena_methods
setGeneric("clusters", function(object, method) standardGeneric("clusters"))


#' @aliases clusters,clena_methods
setMethod("clusters",signature(object="clena"),
          function(object, method=clusterMethods(object)) {
              method <- match.arg(method, clusterMethods(object))
              return(object@clusterObjs[[method]])})


#' mat: get the original data from a clena object.
#' @inheritParams clusterMethods
#' @examples
#' mat(clena_result)
#' @return mat: a matrix
#' @export mat
#' @docType methods
#' @rdname clena_methods
setGeneric("mat", function(object) standardGeneric("mat"))


#' @aliases mat,clena_methods
setMethod("mat",signature(object="clena"),
          function(object) return(object@mat))


#' summary: a summary of a clena object.
#' @return summary: a summary of a clena object.
#' @examples
#' summary(clena_result)
#' @rdname clena_methods
#' @exportMethod summary
setMethod("summary","clena",
          function(object, digits = max(3,getOption("digits")-3)) {
              cat("\nClustering Methods:\n",clusterMethods(object),"\n\n")
              cat("The Number of Clusters:\n",nClusters(object),"\n\n")
              cat("Metric of Distance Matrix:\n", object@metric, "\n\n")
              cat("Agglomeration method for hierarchical clustering (hclust and agnes):\n", object@method, "\n\n")
          })