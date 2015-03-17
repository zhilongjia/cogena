#' get the best clustering methods and the number of clusters
#' 
#' get the best clustering methods and the number of clusters, based that the number
#' of gene sets which are signifigant should be maximum.
#'
#' @inheritParams clusterMethods
#' @param based counting method. Default is "inTotal" to count all the clusters
#' and I, II, All. Other options are "All", "I", "II".
#' @export
#' @docType methods
#' @rdname optCluster
#' @examples
#' data(PD)
#' summary(cogena_result)
#' \dontrun{
#' score <- optCluster(cogena_result)
#' score <- optCluster(cogena_result, based="All")
#' }
#' 
setGeneric("optCluster",  function(object, based="inTotal") standardGeneric("optCluster"))

#' @rdname optCluster
#' @aliases optCluster,cogena
setMethod("optCluster", signature(object="cogena"),
    function(object, based="inTotal"){
    
    based <- match.arg(based, c("inTotal","I", "II", "All"))
    score <- matrix(NA, nrow=length(clusterMethods(object)), ncol=length(nClusters(object)), dimnames=list(clusterMethods(object), nClusters(object)))
    for (i in clusterMethods(object)){
        for (j in as.character(nClusters(object))) {
            
            enrichment_score <- enrichment(object, i, j, roundvalue=FALSE)

            if (is.logical(enrichment_score)) {
                score[i,j]=NA
            } else {
                if (is.logical(object@measures[[i]][[j]])){
                    score[i,j] <- NA
                } else {
                    up_dn_score <- 0
                    #up_dn2 <- 0
                    for (k in 1:j){

                        #print (paste(i, j, k))
                        data <- mat(object)[geneInCluster(object, i, j, as.character(k)),, drop=FALSE]
                        #print (dim(data))

                        up_dn <- apply(data, 1, upORdn, object@sampleLabel)
                        if (length(table(up_dn)) == 1){
                            up_dn_score <- up_dn_score + 1
                        }
                    }

                    up_dn_score <- ifelse (up_dn_score >= as.numeric(j) * 0.75, 1, -1)
                    if (is.null(based)){
                        score[i,j] <- ncol(enrichment_score) * up_dn_score
                    } else if (any(grepl(paste0(based, "#"), rownames(enrichment_score)))) {
                        score[i,j] <- length(enrichment_score[rownames(enrichment_score)[grep(paste0(based, "#"), rownames(enrichment_score))],]) * up_dn_score
                    }

                    up_dn_score <- 0
                }
            }
        }
        # if the 2 clusters have a negative value, all other clusters in this method will have a negative value or NA.
        if (!is.na(score[i,"2"]) && score[i,"2"] < 0){
            for (j in as.character(nClusters(object))){
                if (!is.na(score[i,j]) && score[i,j] >0){
                    score[i,j] <- -score[i,j]
                }
            }
        }
    }
    return (score)
})


upORdn <- function (dat, Label){
    ifelse (mean(dat[which(Label==names(table(Label)[1]))])< mean(dat[which(Label==names(table(Label)[2]))]), 1,-1)
}