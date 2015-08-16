#' heatmap of the gene set enrichment from a cogena object.
#'
#' heatmap of the gene set enrichment score. After obtaining the ennrichemt of 
#' clusters for each gene set, the heatmapPEI will show it as a heatmap in 
#' order. The value shown in heatmapPEI is the -log2(fdr), representing the
#' enrichment score. 
#' 
#' The x-axis shows cluster i and the number of genes in cluster, with
#' red means cluster containing up-regulated genes, green means down-regulated
#' genes, black means there are up and 
#' down regulated genes in this cluster and blue means all DEGs. If parameter
#' add2 is true, another two columns will be shown as well, representing the
#' up and down regulated genes. 
#' 
#' The direction of DEGs are based on latter Vs
#' former from sample labels. For example, labels are 
#' as.factor(c("ct", "Disease")), the "Disease" are latter compared with "ct".
#' Usually, the order is the alphabet.
#' 
#' The y-axis represents the gene sets enriched.
#'
#' @inheritParams enrichment
#' @param low colour for low end of gradient.
#' @param high colour for high end of gradient.
#' @param na.value Colour to use for missing values.
#' @param maintitle a character. like GSExxx. the output of figure will like
#' "cogena: kmeans 3 GSExxx" in two lines. Default is NULL
#' @param printGS print the enriched gene set names or not. Default is TRUE.
#' @param add2 enrichment score for add Up and Down reuglated genes.
#' 
#' @return a gene set enrichment heatmap
#' 
#' @seealso \code{\link{clEnrich}} and \code{\link{heatmapCluster}}
#' 
#' @details
#' orderMethod:
#' \itemize{
#' \item max. ordered by the max value in clusters beside all
#' \item mean. ordered by the mean value in clusters beside all
#' \item All. ordered by all genes
#' \item Up. ordered by up-regulated genes (add2 should be TRUE)
#' \item Down. ordered by down-regulated genes (add2 should be TRUE)
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
#' heatmapPEI(clen_res, "kmeans", "2")
#' heatmapPEI(clen_res, "kmeans", "2", orderMethod="mean")
#' heatmapPEI(clen_res, "kmeans", "3", CutoffNumGeneset=20, 
#'     low = "#132B43", high = "#56B1F7", na.value = "grey50")
#' }
#' 
setGeneric("heatmapPEI", 
    function(object, method, nCluster, CutoffNumGeneset=20,
        CutoffPVal=0.05, orderMethod="max", roundvalue=TRUE,
        low="green", high="red", na.value="white", 
        maintitle=NULL, printGS=TRUE, add2=TRUE)
    standardGeneric("heatmapPEI"))

#' @rdname heatmapPEI
#' @aliases heatmapPEI,cogena
setMethod("heatmapPEI", signature(object="cogena"),
    function(object, method=clusterMethods(object), 
        nCluster=nClusters(object), 
        CutoffNumGeneset=20, CutoffPVal=0.05,
        orderMethod="max", roundvalue=TRUE,
        low="grey", high="red", na.value="white", 
        maintitle=NULL, printGS=TRUE, add2=TRUE) {
        method <- match.arg(method, clusterMethods(object))
        nCluster <- match.arg(nCluster, as.character(nClusters(object)))
        
        enrichment <- enrichment(object, method, nCluster, CutoffNumGeneset, 
            CutoffPVal, orderMethod, roundvalue, add2 =add2)
        
        if (length(enrichment)==1 && is.na(enrichment)){
            return(paste("No enrichment above the cutoff for", method, 
                "when the number of clusters is", nCluster, 
                "with the orderMethod:", orderMethod, "!"))
        }
        if (printGS==TRUE) {
            cat (rev(colnames(enrichment)), sep ="\t")
        }

        cl_color0 <- upDownGene(object, method, nCluster, add2)
        
        cl_color <- sapply(cl_color0, function(x) {if (x>=0.6) "red" else if (x<=-0.6) "green4" else "black"})
        cl_color <- c(cl_color, "blue")
        

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

