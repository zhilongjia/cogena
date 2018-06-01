#' paralleled clustering
#'
#' Internal functions.
#'
#' Paralleled clustering with different clustering methods and different number
#' of clusters. These functions should not be called directly by the user.
#' @return a list containing the cluster information 
#' @import amap
#' @import cluster
#' @importFrom mclust Mclust mclustBIC
#' @import fastcluster
#' @import parallel
#' @importFrom kohonen som somgrid
#' @import foreach
#' @import doParallel
#' @keywords internal
#' 
pClusters <- function(mat, Distmat, clMethod, nClust, method, 
    metric, ncore, verbose, ...) {

    #cluster init
    switch(clMethod,
        hierarchical = { clusterObj <- fastcluster::hclust(Distmat,method)},
        diana = { clusterObj <- diana(Distmat, ...)},
        kmeans = {
            clusterObj <- vector("list",length=length(nClust))
            names(clusterObj) <- nClust
            clusterObjInit <- fastcluster::hclust(Distmat, method)},
        agnes = {clusterObj <- agnes(Distmat, method=method, ...)},
        apcluster = {
            ap_cluster_infor <- vector("integer", length=length(rownames(mat)))
            names(ap_cluster_infor) <- rownames(mat)
            clusterObj0 <- apcluster::apcluster(apcluster::negDistMat(r=2), mat)
            clusterObj <- list()
        },
        {clusterObj <- list()})

    # parallel the Clustering with the number of Cluster is nc
    #doMC::registerDoMC(ncore)
    cl <- parallel::makeCluster(ncore)
    doParallel::registerDoParallel(cl)
    #doParallel::registerDoParallel(cores=ncore)
    
    if (verbose) {print(paste("getDoParWorkers:", foreach::getDoParWorkers()))}
    nc=NULL
    nClust <- as.numeric(nClust)
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
                clusterObj <- amap::Kmeans(mat, centers=initial, 
                    iter.max=100, 
                    method=ifelse((metric=="NMI" || metric=="biwt"), 
                    "correlation", metric),...)
            },
            fanny = {
                clusterObj <- cluster::fanny(Distmat, nc, ...)
                clusterObj$cluster <- clusterObj$clustering
            },
            model = {
                clusterObj <- mclust::Mclust(mat,nc, ...)
                clusterObj$cluster <- clusterObj$classification
            },
            som = {
                clusterObj <- kohonen::som(mat, grid=kohonen::somgrid(1,nc), ...)
                clusterObj$cluster <- clusterObj$unit.classif
                names(clusterObj$cluster) <- rownames(mat)
            },
            pam = {
                clusterObj <- pam(Distmat, nc, ...)
                clusterObj$cluster <- clusterObj$clustering
            },
            clara = {
                clusterObj <- clara(mat, nc, 
                    metric=ifelse(metric=="manhattan", metric, "euclidean"),
                    ...)
                clusterObj$cluster <- clusterObj$clustering
            },
            sota = {
                suppressWarnings( clusterObj <- sota(mat, nc-1, 
                    distance=ifelse(metric=="euclidean", metric, 
                        "correlation")) )
            },
            apcluster = {
                if (nc <= length(clusterObj0) ) {
                    ap_cluster <- cutree(apcluster::as.hclust(clusterObj0), nc)
                    for (i in 1:nc) {
                        ap_cluster_infor[names(unlist(clusterObj0@clusters[which(ap_cluster == i)]))] = i
                    }
                    clusterObj$cluster <- ap_cluster_infor
                } else {
                    clusterObj$cluster <- NA
                }
                
            },
            ## otherwise - hierarchical, diana, agnes
            {clusterObj$cluster <- cutree(as.hclust(clusterObj), nc)}
            )

    clusterObj$cluster
    } #END OF NC LOOP
    
    # stopImplicitCluster()
    parallel::stopCluster(cl)
    
    
    return (clusterList)
}

