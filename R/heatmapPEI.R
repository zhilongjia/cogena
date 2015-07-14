#' heatmap of the gene set enrichment from a cogena object.
#'
#' heatmap of the gene set enrichment index. After obtaining the ennrichemt of 
#' clusters in the gene sets, the heatmapPEI will show it as a heatmap with 
#' order.
#'
#' @inheritParams enrichment
#' @param low colour for low end of gradient.
#' @param high colour for high end of gradient.
#' @param na.value Colour to use for missing values.
#' @param maintitle a character. like GSExxx. the output of figure will like
#' "cogena: kmeans 3 GSExxx" in two lines. Default is NULL
#' @param printGS print the enriched gene set names or not. Default is TRUE.
#' @return a gene set enrichment heatmap
#' @seealso \code{\link{clEnrich}} and \code{\link{heatmapCluster}}
#' 
#' @details
#' orderMethod:
#' \itemize{
#' \item max. ordered by the max value in clusters beside all
#' \item mean. ordered by the mean value in clusters beside all
#' \item All. ordered by all genes
#' \item I. ordered by the I cluster in only two clusters 
#' (Up or Down-regulated)
#' \item II. ordered by the II cluster in only two clusters 
#' (Up or Down-regulated)
#' }
#' 
#' @export
#' @import ggplot2
#' @import reshape2
#' @docType methods
#' @rdname heatmapPEI
#' @examples
#' data(PD)
#' annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt", 
#' package="cogena")
#' 
#' \dontrun{
#' genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=c("hierarchical","kmeans"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
#' 
#' #summay this cogena object
#' summary(clen_res)
#' 
#' #heatmapPEI
#' heatmapPEI(clen_res, "kmeans", "2", orderMethod="mean")
#' heatmapPEI(clen_res, "kmeans", "3", CutoffNumGeneset=20, 
#'     low = "#132B43", high = "#56B1F7", na.value = "grey50")
#' }
#' 
setGeneric("heatmapPEI", 
    function(object, method, nCluster, CutoffNumGeneset=20,
        CutoffPVal=0.05, orderMethod="max", roundvalue=TRUE,
        low="green", high="red", na.value="white", 
        maintitle=NULL, printGS=TRUE) 
    standardGeneric("heatmapPEI"))

#' @rdname heatmapPEI
#' @aliases heatmapPEI,cogena
setMethod("heatmapPEI", signature(object="cogena"),
    function(object, method=clusterMethods(object), 
        nCluster=nClusters(object), 
        CutoffNumGeneset=20, CutoffPVal=0.05,
        orderMethod="max", roundvalue=TRUE,
        low="grey", high="red", na.value="white", 
        maintitle=NULL, printGS=TRUE) {
        method <- match.arg(method, clusterMethods(object))
        nCluster <- match.arg(nCluster, as.character(nClusters(object)))
        
        enrichment <- enrichment(object, method, nCluster, CutoffNumGeneset, 
            CutoffPVal, orderMethod, roundvalue)
        if (length(enrichment)==1 && is.na(enrichment)){
            return(paste("No enrichment above the cutoff for", method, 
                "when the number of clusters is", nCluster, "!"))
        }
        if (printGS==TRUE) {
            cat (rev(colnames(enrichment)), sep ="\t")
        }

        if (nCluster != "2"){
            cl_color <- c(rep("black", nrow(enrichment)-3), rep("blue", 3))
        } else {
            cl_color <- rep("black", 3)
        }

        enrichment <- reshape2::melt(enrichment)
        #legend breaks
        if (max(enrichment$value, na.rm=TRUE) > 15 ){
            breaks <- seq(15, max(enrichment$value, na.rm=TRUE), 10)
        } else {
            breaks <- NULL
        }
        Var1=Var2=value=NULL
        if (!is.null(title)) {
            title=paste("cogena:", method, nCluster, "\n", maintitle)
        } else {
            title=paste("cogena:", method, nCluster)
        }

        ggplot2::ggplot(enrichment, aes(as.factor(Var1), Var2)) + 
            geom_tile(aes(fill = value)) + 
            scale_fill_gradient2("score",  mid=low, midpoint=4, low=low, 
                high=high, na.value=na.value, breaks=c(4.32, breaks)) +
            geom_text(aes(fill=value, label=value),size=4, na.rm=TRUE) +
            labs(list(title = title, x = "Cluster", y = "Gene set")) +
            theme(axis.text.y = element_text(size = rel(1.5), face="bold")) +
            theme(axis.text.x = element_text(size = rel(1.3), angle=30, 
                face="bold", color=cl_color))
    }
)

