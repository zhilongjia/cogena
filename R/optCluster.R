#' get the best clustering methods and the number of clusters
#' 
#' get the best clustering methods and the number of clusters, based that 
#' the number of gene sets which are signifigant should be maximum.
#'
#' @inheritParams clusterMethods
#' @param based counting method. Default is "inTotal" to count all enriched 
#' gene sets from all the clusters and I, II, All. Other options are "All",
#'  "I", "II". Here "I" or "II" means all the up-regualted or down-regualted
#'  genes (Do check the output of function heatmapCluster to suggest which is
#'  up-regulated (down-regulated)!), while "All" means all the DEGs used.
#' @param ncores cores used for caculating optCluster. Default is same as ncores
#' used during cogena function, but it will be the same as number of cores 
#' machine has if ncores parameter is exceed it.
#' @param CutoffPVal the cut-off of p-value. Default is 0.05.
#' @return a score matrix
#' @export
#' @docType methods
#' @rdname optCluster
#' @import parallel
#' @import foreach
#' @import doParallel
#' @examples
#' data(PD)
#' annofile <- system.file("extdata", "c2.cp.kegg.v4.0.symbols.gmt", 
#' package="cogena")
#' cogena_result <- cogena(DEexprs, nClust=2:3, 
#' clMethods=c("hierarchical","kmeans"), metric="correlation", 
#' method="complete",  annofile=annofile, sampleLabel=sampleLabel, 
#' ncore=1, verbose=TRUE)
#' summary(cogena_result)
#' 
#' score <- optCluster(cogena_result)
#' score <- optCluster(cogena_result, based="All")
#' 
#' 
setGeneric("optCluster", 
    function(object, based="inTotal", ncores=object@ncore,
        CutoffPVal=0.05) standardGeneric("optCluster"))

#' @rdname optCluster
#' @aliases optCluster,cogena
setMethod("optCluster", signature(object="cogena"),
    function(object, based="inTotal", ncores=object@ncore, CutoffPVal=0.05){

    based <- match.arg(based, c("inTotal","I", "II", "All"))

    if (ncores > parallel::detectCores()) {
        ncores <- parallel::detectCores()
    }
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    # doParallel::registerDoParallel(cores=ncores)
    i = NULL
    score <- 
        foreach::foreach(i = clusterMethods(object), .combine='rbind') %dopar% {
        cl_score <- vector(mode="numeric", length = length(nClusters(object)))
        names(cl_score) <- as.character(nClusters(object))
        cl_score_index = 1
        for (j in as.character(nClusters(object))) {
            enrichment_score <- enrichment(object, i, j, roundvalue=FALSE)

            if (is.logical(enrichment_score) || 
                is.logical(object@measures[[i]][[j]])) {
                cl_score[cl_score_index] = NA
            } else {
                up_dn_score <- 0
                for (k in 1:j){
                    data <- mat(object)[geneInCluster(object, i, j, 
                        as.character(k)),, drop=FALSE]
                    up_dn <- apply(data, 1, upORdn, object@sampleLabel)
                    if (length(table(up_dn)) == 1){
                        up_dn_score <- up_dn_score + 1
                    }
                }

                up_dn_score <- 
                    ifelse (up_dn_score >= as.numeric(j) * 0.75, 1, -1)
                if (based=="inTotal"){
                    cl_score[cl_score_index] <- 
                        ncol(enrichment_score) * up_dn_score
                } else if (
                    any(grepl(paste0(based, "#"), rownames(enrichment_score)))){
                    cl_score[cl_score_index] <- 
                        length(which(enrichment_score[rownames(enrichment_score)
                            [grep(paste0("^", based, "#"), 
                                rownames(enrichment_score))],] > 
                                    -log2(CutoffPVal))) * up_dn_score
                }
                up_dn_score <- 0
            }
            cl_score_index = cl_score_index + 1
        }
        
        # if the 2 clusters have a negative value, all other clusters 
        #in this method will have a negative value or NA.
        if (!is.na(cl_score["2"]) && cl_score["2"] < 0) {
            for (j in as.character(nClusters(object))){
                if (!is.na(cl_score[j]) && cl_score[j] >0){
                    cl_score[j] <- -cl_score[j]
                }
            }
        }
        cl_score
        }
    stopCluster(cl)
    # stopImplicitCluster()
    
    rownames(score) <- clusterMethods(object)
    return (score)
})


upORdn <- function (dat, Label){
    ifelse (mean(dat[which(Label==names(table(Label)[1]))])< 
        mean(dat[which(Label==names(table(Label)[2]))]), 1,-1)
}
