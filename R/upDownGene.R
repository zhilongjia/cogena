#' Show up or down-regualted genes for a clustering method and the number of clusters.
#'
#' The value means up or down regulated genes for each cluster. 1 suggests that
#' genes in the cluster is up-regualted genes, while -1 down-regualted genes. 
#' value within (-1, 1) means genes there are both up and down regulated genes
#' in the cluster. Return a vector with the length of nCluster if add2 is 
#' FALSE, or the length of nCluster + 2 if add2 is TRUE and nCluster is not 2. 
#' In the latter situation, the last two itemes represent Up and Down
#' reuglated genes
#' 
#' @param object a genecl or cogena object
#' @param method as clMethods in genecl function
#' @param nCluster cluster number
#' @param add2 add2 enrichment score for add Up and Down reuglated genes
#' @return a vector
#' @rdname upDownGene
#' @export
#' 
#' @examples 
#' data(PD)
#' annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt", 
#' package="cogena")
#' 
#' 
#' genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
#' 
#' upDownGene(clen_res, "kmeans", "3", add2=TRUE)
#' 
#' upDownGene(clen_res, "kmeans", "2", add2=FALSE)
#' 
setGeneric("upDownGene", 
           function(object, method, nCluster, add2=FALSE) 
               standardGeneric("upDownGene"))

#' @rdname upDownGene
#' @aliases upDownGene,cogena
setMethod("upDownGene", signature(object="cogena"), function(
    object, method, nCluster, add2=FALSE) {
    sampleLabel <- object@sampleLabel
    
    logFC <- NULL; cluster_id = NULL
    
    # Return a vector representing the up or down regulated genes. (1 or -1)
    logfc <- function(object, method, nClusterX, sampleLabel) {
        geneExp <- as.data.frame(geneExpInCluster(object, method, nClusterX)$clusterGeneExp)
        sampleLabel <- geneExpInCluster(object, method, nClusterX)$label
        
        geneExp$meanA <- apply(geneExp[,names(sampleLabel)[which((sampleLabel == names(table(sampleLabel))[1]))]], 1, mean)
        geneExp$meanB <- apply(geneExp[,names(sampleLabel)[which((sampleLabel == names(table(sampleLabel))[2]))]], 1, mean)
        geneExp$logFC <- ifelse( log2(geneExp$meanB/geneExp$meanA)>0, 1, -1)
        
        geneExp <- geneExp[,c("cluster_id", "logFC")]
        cluster_upDn <- dplyr::group_by(geneExp, cluster_id)
        cluster_upDn <- dplyr::summarise(cluster_upDn, UpDn=sum(logFC)/n() )
        cluster_upDn <- dplyr::arrange(cluster_upDn, cluster_id)
        cluster_upDn$UpDn
    }
    cluster_logfc <- logfc(object, method, nCluster)
    
    if (isTRUE(add2) ) {
        res <- c(cluster_logfc, c(1, -1))
    } else {
        res <- cluster_logfc
    }
    return (res)
})
    

