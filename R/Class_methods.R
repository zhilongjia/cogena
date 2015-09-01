#' Basic methods for a genecl object.
#' 
#' clusterMethods: get the methods of clustering used.
#' @param object a genecl or cogena object
#' @exportMethod clusterMethods
#' @import methods
#' @docType methods
#' @rdname clusterMethods
#' @return clusterMethods: a character vector.
#' @examples
#' data(Psoriasis)
#' genecl_result <- coExp(DEexprs, nClust=2:3, 
#' clMethods=c("hierarchical","kmeans"), metric="correlation", 
#' method="complete", 
#' ncore=1, verbose=TRUE)
#' clusterMethods(genecl_result)
#' 
#' 
setGeneric("clusterMethods", 
    function(object) standardGeneric("clusterMethods"))

#' @rdname clusterMethods
#' @aliases clusterMethods,genecl_methods
setMethod("clusterMethods",signature(object="genecl"),
    function(object) return(object@clMethods))

#' @rdname clusterMethods
#' @aliases clusterMethods,cogena_methods
setMethod("clusterMethods",signature(object="cogena"),
          function(object) return(object@clMethods))

################################################################################
#' nClusters: get the number of clusters from a genecl object.
#' @inheritParams clusterMethods
#' @examples
#' \dontrun{
#' nClusters(genecl_result)
#' }
#' @return nClusters: a numeric vector.
#' @exportMethod nClusters
#' @docType methods
#' @rdname nClusters
setGeneric("nClusters", function(object) standardGeneric("nClusters"))

#' @rdname nClusters
#' @aliases nClusters,genecl_methods
setMethod("nClusters",signature(object="genecl"),
    function(object) return(object@nClust))

#' @rdname nClusters
#' @aliases nClusters,cogena_methods
setMethod("nClusters",signature(object="cogena"),
          function(object) return(object@nClust))

################################################################################
#' geneclusters: get the cluster information of a certain clustering method with
#' a certain number.
#' @inheritParams clusterMethods
#' @param method as clMethods in genecl function
#' @param nClust cluster numbers
#' @return geneclusters: a list or hclust depends on the method
#' @examples
#' \dontrun{
#' geneclusters(genecl_result, "kmeans", 3)
#' geneclusters(genecl_result, "hierarchical", 4)
#' }
#' @exportMethod geneclusters
#' @docType methods
#' @rdname geneclusters
setGeneric("geneclusters", function(object, method, nClust) standardGeneric("geneclusters"))

#' @rdname geneclusters
#' @aliases geneclusters,genecl_methods
setMethod("geneclusters",signature(object="genecl"),
    function(object, method=clusterMethods(object), nClust) {
        method <- match.arg(method, clusterMethods(object))
        nClust <- as.character(nClust)
        return(object@clusterObjs[[method]][[nClust]]) })

#' @rdname geneclusters
#' @aliases geneclusters,cogena_methods
setMethod("geneclusters",signature(object="cogena"),
          function(object, method=clusterMethods(object), nClust) {
              method <- match.arg(method, clusterMethods(object))
              nClust <- as.character(nClust)
              return(object@clusterObjs[[method]][[nClust]]) })

################################################################################
#' mat: get the original data from a genecl object.
#' @inheritParams clusterMethods
#' @examples
#' \dontrun{
#' mat(genecl_result)
#' }
#' @return mat: a matrix
#' @exportMethod mat
#' @docType methods
#' @rdname mat
setGeneric("mat", function(object) standardGeneric("mat"))

#' @rdname mat
#' @aliases mat,genecl_methods
setMethod("mat",signature(object="genecl"),
    function(object) return(object@mat))

#' @rdname mat
#' @aliases mat,cogena_methods
setMethod("mat",signature(object="cogena"),
          function(object) return(object@mat))


################################################################################
#' summary: a summary of a genecl object.
#' @inheritParams clusterMethods
#' @return summary: a summary of a genecl object.
#' @examples
#' \dontrun{
#' summary(genecl_result)
#' }
#' @rdname summary
#' @exportMethod summary
setMethod("summary","genecl",
    function(object) {
    cat("\nClustering Methods:\n",clusterMethods(object),"\n\n")
    cat("The Number of Clusters:\n",nClusters(object),"\n\n")
    cat("Metric of Distance Matrix:\n", object@metric, "\n\n")
    cat("Agglomeration method for hierarchical 
        clustering (hclust and agnes):\n", object@method, "\n\n")
    }
)

#' @rdname summary
#' @aliases summary,cogena_methods
setMethod("summary","cogena",
          function(object) {
              cat("\nClustering Methods:\n",clusterMethods(object),"\n\n")
              cat("The Number of Clusters:\n",nClusters(object),"\n\n")
              cat("Metric of Distance Matrix:\n", object@metric, "\n\n")
              cat("Agglomeration method for hierarchical clustering (hclust and agnes):\n", object@method, "\n\n")
              cat("Gene set:\n", object@gmt, "\n\n")
          }
)

################################################################################
#' show: show the class of cogena or genecl object
#' @param object a genecl or cogena object
#' @examples
#' \dontrun{
#' show(genecl_result)
#' }
#' @rdname show
#' @exportMethod show
setMethod(f="show", signature="cogena",
    function(object) { 
        cat("An instance of ", "\"", class(object), "\".", "\n", sep="")
})

#' @rdname show
#' @aliases show,cogena_methods
#' @return show which instance
setMethod(f="show", signature="genecl",
    function(object) { 
        cat("An instance of ", "\"", class(object), "\".", "\n", sep="")
})

