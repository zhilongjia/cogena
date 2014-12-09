#' get the best clustering methods and the number of clusters
#' 
#' get the best clustering methods and the number of clusters, based that the number
#' of gene sets which are signifigant should be maxium.
#'
#' @inheritParams clusterMethods
#' @export
#' @docType methods
#' @rdname optClustering
#' @examples
#' ##
setGeneric("optClustering", function(object) standardGeneric("optClustering"))


#' @aliases optClustering, clena
setMethod("optClustering", signature(object="clena"),
    function(object){
    score <- matrix(NA, nrow=length(clusterMethods(object)), ncol=length(nClusters(object)), dimnames=list(clusterMethods(object), nClusters(object)))
    for (i in clusterMethods(object)){
        for (j in as.character(nClusters(object))) {
            enrichment_score <- enrichment(object, i, j, roundvalue=FALSE)
            if (is.logical(enrichment_score)) {
                score[i,j]=NA
            } else {
                score[i,j] <- ncol(enrichment_score)
            }
        }
    }
    return (score)
})


