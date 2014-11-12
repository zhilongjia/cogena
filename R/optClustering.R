#! /usr/bin/Rscript

# require(devtools)
# load_all("~/clena", export_all=FALSE)
# 
# object <-  GSE20292.clena.cluster

optClustering <- function(object){
    score <- matrix(NA, nrow=length(clusterMethods(object)), ncol=length(nClusters(object)), dimnames=list(clusterMethods(object), nClusters(object)))
    for (i in clusterMethods(object)){
        for (j in as.character(nClusters(object))) {
            enrichment_score <- enrichment0(object, i, j, roundvalue=FALSE)
            if (is.logical(enrichment_score)) {
                score[i,j]=NA
            } else {
                score[i,j] <- length(which(apply(enrichment_score, 2, max, na.rm=TRUE)>=-log2(0.05)))
            }
        }
    }
    return (score)
}


#tmp <- optClustering(GSE20292.clena.cluster)
# i = clusterMethods(object)[1]
# j = nClusters(object)[1]

