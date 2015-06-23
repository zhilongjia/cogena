#' heatmap of gene expression profilings with cluster indication.
#' 
#' heatmap of gene expression profilings with cluster-based color indication.
#' @inheritParams enrichment
#' @param sampleColor a color vector with the sample length. 
#' The default is c("darkblue", "cyan"). Or setting as "random" if 
#'  sample Label has more than 2 levels.
#' @param clusterColor a color vector with the cluster length. 
#' The default is rainbow(nClusters(object)).
#' @param clusterColor2 a color vector with 2 elements. The default is 
#' c("coral3", "deepskyblue1").
#' @param heatmapcol col for heatmap. The default is greenred(75).
#' @param maintitle a character. like GSExxx. the output of figure will like
#' "kmeans 3 Clusters GSExxx" in two lines.
#' @param printSum print the summary of the number of genes in each cluster. 
#' Default is TRUE.
#' @param ... other parameters to heatmap.3.
#' @return a gene expression heatmap with Cluster information figure
#' @export
#' @import gplots
#' @rdname heatmapCluster
#' @docType methods
#' @seealso \code{\link{cogena}}, \code{\link{heatmap.3}} and
#' \code{\link{heatmapPEI}}
#' @examples
#' data(PD)
#' annofile <- system.file("extdata", "c2.cp.kegg.v4.0.symbols.gmt", 
#' package="cogena")
#' cogena_result <- cogena(DEexprs, nClust=2:3, 
#' clMethods=c("hierarchical","kmeans"), metric="correlation", 
#' method="complete",  annofile=annofile, sampleLabel=sampleLabel, 
#' ncore=1, verbose=TRUE)
#' 
#' #summay this cogena object
#' summary(cogena_result)
#'
#' #heatmapCluster
#' 
#' heatmapCluster(cogena_result, "hierarchical", "3")
#' heatmapcol <- gplots::redgreen(75) 
#' heatmapCluster(cogena_result, "hierarchical", "3", heatmapcol=heatmapcol)
#' 
#' # Setting sampleColor as "random" if sample Label has more 
#' # than 2 levels.
#' heatmapCluster(cogena_result, "hierarchical", "3", sampleColor="random")
#' 
setGeneric("heatmapCluster", 
    function(object, method, nClusters, sampleColor=c("darkblue", "cyan"),
        clusterColor=NULL, clusterColor2=NULL, heatmapcol=NULL, maintitle=NULL,
        printSum=TRUE, ...) 
    standardGeneric("heatmapCluster"))

#' @rdname heatmapCluster
#' @aliases heatmapCluster
setMethod("heatmapCluster", signature(object="cogena"),
    function (object, method=clusterMethods(object), 
        nClusters=nClusters(object),
        sampleColor=c("darkblue", "cyan"), clusterColor=NULL,
        clusterColor2=NULL, heatmapcol=NULL, maintitle=NULL, 
        printSum=TRUE, ...){

    #get the parameters
    method <- match.arg(method, clusterMethods(object))
    nClusters <- match.arg(nClusters, as.character(nClusters(object)))
    mat <- mat(object)

    if (method %in% c("hierarchical", "diana", "agnes")) {
        cluster_size <- cutree(clusters(object, method), nClusters)
        names(cluster_size) <- rownames(mat)
        cluster_size2 <- cutree(clusters(object, method), "2")
        names(cluster_size2) <- rownames(mat)
    } else {
        cluster_size <- clusters(object, method)[[nClusters]]$cluster
        cluster_size2 <- clusters(object, method)[["2"]]$cluster
    }

    if (printSum==TRUE) {
        cat ("The number of genes in each cluster:\n")
        print (table(cluster_size2))
        if (nClusters != "2"){
            print (table(cluster_size))
        }
    }

    #print (cluster_size)
    #reorder the mat based on the clustering and type of sample.
    mat <- mat[order(cluster_size, decreasing=FALSE), order(object@sampleLabel)]
    #add the type of sample into the colnames
    #colnames(mat) <- paste(colnames(mat), sampleLabel, sep="_")
    
    #color setting
    sampleLabel <- sort(object@sampleLabel)
    # ColSideColors <- map2col(as.numeric(as.factor(sampleLabel)), sampleColor)
    if (sampleColor == "random") {
        sampleColor <- sample(topo.colors(nlevels(as.factor(sampleLabel))))
    }
    ColSideColors <- map2col(as.numeric(as.factor(sampleLabel)), sampleColor)
    
    if (is.null(clusterColor)) {
        clusterColor <- sample(rainbow(nClusters)) #, alpha = c(1, 0.6)
    }
    if (is.null(clusterColor2)){
        clusterColor2 <- c("coral3", "deepskyblue1")
    }

    #clusterColor <- rainbow(nClusters, alpha = c(1, 0.6))
    #make the neighboor colors different
    #clusterColor <- clusterColor[clusterColor]
    
    RowSideColors <- map2col(sort(cluster_size, decreasing=FALSE), clusterColor)
    RowSideColors2 <- map2col(cluster_size2[names(sort(cluster_size, 
        decreasing=FALSE))], clusterColor2)
    if (is.null(heatmapcol)) {
        heatmapcol <- greenred(75)
    }
    if (nClusters != "2"){
        RowSideColors <- t(cbind(RowSideColors2, RowSideColors))
        rownames(RowSideColors) <- paste("Size:", c(2, nClusters))
    } else {
        RowSideColors <- t(as.matrix(RowSideColors))
    }

    if (!is.null(maintitle)) {
        maintitle=paste(method, nClusters, "clusters", "\n", maintitle)
    } else {
        maintitle=paste(method, nClusters, "clusters",sep="_")
    }
    
    heatmap.3(mat, col=heatmapcol, trace="none", scale="row", Rowv=FALSE, 
        Colv=FALSE, dendrogram="none", labRow=NA, 
        # colsep=length(which(sampleLabel==sampleLabel[1])), 
        colsep=cumsum( table(sampleLabel)),
        rowsep=cumsum( table(cluster_size)), #adjCol=c(0.8,0),
        sepcolor="white", sepwidth=c(0.05,1),
        key=TRUE, symkey=FALSE, density.info="none", keysize=1.5, 
        main=maintitle,
        ColSideColors=as.matrix(ColSideColors), ylab="Genes", cexCol=0.8,
        RowSideColors=RowSideColors, 
        RowSideColorsSize=ifelse(nClusters!=2,2,1),
        lhei=c(1.2, 4), ...)
        par(lend = 0, xpd=TRUE)
        legend("left", legend = paste0(1:nClusters),
            col = clusterColor, lty= 1, lwd = 20, bty = "n", 
            title = paste(nClusters, "Clusters"))
        if (nClusters != "2"){
            legend("bottomleft", legend = as.character(as.roman(1:2)),
                col = clusterColor2, lty= 1, lwd = 20, bty = "n", 
                title = paste(2, "Clusters"))
        }
        legend("top", legend = names(table(sampleLabel)), col = sampleColor, 
            lty=1, lwd=20, bty = "n", title = "Type of Sample", horiz=TRUE)
    # return (mat)
    })


################################################################################
#map2color: get the color vector from the numeric vector x using pal, 
#such as, rainbow(200)
#Example: 
#map2color(cutree(clusters(cogena.cluster, "hierarchical"), 3),rainbow(200))
#http://stackoverflow.com/questions/15006211/
#how-do-i-generate-a-mapping-from-numbers-to-colors-in-r
##############################################################################
map2col<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), 
        all.inside=TRUE)]
}
