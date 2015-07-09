#' paralleled clustering
#'
#' Internal cogena functions.
#'
#' Paralleled clustering with different clustering methods and different number
#' of clusters. These functions should not be called directly by the user.
#' @return a list containing the cluster information 
#' @import amap
#' @import cluster
#' @import mclust
#' @import fastcluster
#' @import parallel
#' @import foreach
#' @import doParallel
#' @keywords internal
#' 
pClusters <- function(mat, Distmat, clMethod, nClust, method, 
    metric, annotation, ncore, annotationGenesPop, verbose, ...) {

    meausres <- matrix(NA, nrow=length(nClust), ncol=ncol(annotation))


    #cluster init
    switch(clMethod,
        hierarchical = { clusterObj <- fastcluster::hclust(Distmat,method)},
        diana = { clusterObj <- diana(Distmat, ...)},
        kmeans = {
            clusterObj <- vector("list",length=length(nClust))
            names(clusterObj) <- nClust
            clusterObjInit <- fastcluster::hclust(Distmat, method)},
        agnes = {clusterObj <- agnes(Distmat, method=method, ...)},
        {clusterObj <- vector("list",length=length(nClust))
        names(clusterObj) <- nClust })

    # parallel the Clustering with the number of Cluster is nc
    #doMC::registerDoMC(ncore)
    cl <- parallel::makeCluster(ncore)
    doParallel::registerDoParallel(cl)
    #doParallel::registerDoParallel(cores=ncore)
    
    if (verbose) {print(paste("getDoParWorkers:", foreach::getDoParWorkers()))}
    nc=NULL
    clusterList <- foreach::foreach (nc = nClust) %dopar% {
        if (verbose) {print (paste(clMethod, "Starting nClust:", nc))}
        switch(clMethod,
            kmeans = {
                initial <- tapply(mat, 
                    list(rep(cutree(clusterObjInit,nc),ncol(mat)),col(mat)),
                        function(x) mean(x, na.rm=TRUE))
                if(length(dup <- which(duplicated(initial)))>0) {
                    for(dupi in dup) 
                        initial[dupi,] <- 
                            initial[dupi,] + jitter(initial[dupi,])
                }
                dimnames(initial) <- list(NULL,dimnames(mat)[[2]])

                #amap::Kmeans can use other distance metric besides euclidean.
                clusterObj <- Kmeans(mat, centers=initial, 
                    iter.max=100, 
                    method=ifelse((metric=="NMI" || metric=="biwt"), 
                    "correlation", metric),...)
                cluster <- clusterObj$cluster
            },
            fanny = {
                clusterObj <- fanny(Distmat, nc, ...)
                clusterObj$cluster <- clusterObj$clustering
                cluster <- clusterObj$cluster
            },
            model = {
                clusterObj <- mclust::Mclust(mat,nc, ...)
                clusterObj$cluster <- clusterObj$classification
                cluster <- clusterObj$cluster
            },
            som = {
                #              clusterObj <- som(mat, grid=somgrid(1,nc))
                clusterObj <- kohonen::som(mat, grid=somgrid(1,nc), ...)
                clusterObj$cluster <- clusterObj$unit.classif
                names(clusterObj$cluster) <- rownames(mat)
                cluster <- clusterObj$cluster
            },
            pam = {
                #              clusterObj <- pam(Distmat, nc)
                clusterObj <- pam(Distmat, nc, ...)
                clusterObj$cluster <- clusterObj$clustering
                cluster <- clusterObj$cluster
            },
            clara = {
                clusterObj <- clara(mat, nc, 
                    metric=ifelse(metric=="manhattan", metric, "euclidean"),
                    ...)
                clusterObj$cluster <- clusterObj$clustering
                cluster <- clusterObj$cluster
            },
            sota = {
                clusterObj <- sota(mat, nc-1, 
                    distance=ifelse(metric=="euclidean", metric, 
                        "correlation"))
                cluster <- clusterObj$cluster
            },
            ## otherwise - hierarchical, diana, agnes
            {cluster <- cutree(clusterObj, nc)})

        #if (!exists("cluster")) {cluster <- clusterObj$cluster}
        #if (!is.vector(cluster)) {cluster <- clusterObj$cluster}
        if (is.null(names(cluster))) {names(cluster) <- rownames(mat)}

        if (verbose) {print (paste(clMethod, "nClust:", nc, "End"))}

        ## Gene sets enrichment measures
        if (nc != length(unique(cluster))) {
            pei=NA
            warning (paste("Cluster", nc, "(aim) only have", 
                length(unique(cluster)), "(result) clusters"))
            return (list(clusterObj=clusterObj, measures=pei))
        } else {
            pei <- matrix(NA, nrow=length(unique(cluster)), 
                ncol=ncol(annotation))
            rownames(pei) <- sort(unique(cluster))
            colnames(pei) <- colnames(annotation)
        for (k in  sort(unique(cluster))) {
            genenames <- names(which(cluster==k))
            pei[as.character(k),] <- PEI(genenames, annotation=annotation, 
                annotationGenesPop=annotationGenesPop)}
        }

        #calculate the enrichment of the All genes of all cluster
        #if (verbose) print ("pei cluster done")
        All <- PEI(names(cluster), annotation=annotation, 
            annotationGenesPop=annotationGenesPop)
        pei <- rbind(pei, All)

        #negative log p value
        logAdjPEI <- function (pei, verbose=TRUE) {
            #fdr based on pval (pei above)
            pei.adjust <- matrix(p.adjust(pei, "fdr"), ncol=ncol(pei))
            dimnames(pei.adjust) <- dimnames(pei)
            pei.NeglogPval <- -log(pei.adjust)
        }
        pei <- logAdjPEI(pei)
        if(verbose) print(paste("End PEI analysis,", clMethod, nc, "clusters"))
        list(clusterObj=clusterObj, measures=pei)
    } #END OF NC LOOP
    
    stopCluster(cl)
    #stopImplicitCluster()
    
########################################################
    #Combine the paralled results
    measuresComb <- vector("list", length(nClust))
    names(measuresComb) <- nClust

    for (j in seq(length(clusterList))){
        if (clMethod %in% c("hierarchical", "diana", "agnes")) {
            if (j==1){
                clusterObj <- clusterList[[j]]$clusterObj}
        } else {
            clusterObj[[j]] <- clusterList[[j]]$clusterObj
        }
    measuresComb[[j]] <- clusterList[[j]]$measures }

list(clusterObj=clusterObj, measures=measuresComb)
}


