#' heatmap designed for for CMap gene set only
#' 
#' heatmapCmap is desgined for the cogena result from CMap
#' only so as to collapse the multi-isntance drugs in CMap! 
#' 
#' orderMethod:
#' \itemize{
#' \item max. ordered by the max value in clusters beside all
#' \item mean. ordered by the mean value in clusters beside all
#' \item All. ordered by all genes
#' \item Up. ordered by up-regulated genes (add2 should be TRUE)
#' \item Down. ordered by down-regulated genes (add2 should be TRUE)
#' }
#' 
#' MultiInstance:
#' \itemize{
#' \item drug. merge based on cmap_name
#' \item celldrug. merge based on cmap_name and cell type
#' \item conccelldrug. merge based on cmap_name, cell type and concentration
#' }
#' 
#' @inheritParams enrichment
#' @param nCluster as nClust in cogena function.
#' @param orderMethod the order method, default is max, other options are 
#' "mean", "all", "I", "II" or a number meaning the ith cluster.
#' @param MultiInstance merge multi instances. Options are "drug", "celldrug", "conccelldrug", "concdrug".
#' @param CutoffNumGeneset the cut-off of the number of gene sets in the 
#' return table. The default is 20.
#' @param CutoffPVal the cut-off of p-value. The default is 0.05.
#' @param mergeMethod max or mean. The default is mean.
#' @inheritParams heatmapPEI
#' @param maintitle a character. Default is null
#' @param printGS print the enriched gene set names or not. Default is FALSE
#' @param geom tile or circle
#' 
#' @return a gene set enrichment heatmap
#' 
#' 
#' @examples
#' data(Psoriasis)
#' annofile <- system.file("extdata", "CmapDn100.gmt.xz", package="cogena")
#' \dontrun{
#' genecl_result <- coExp(DEexprs, nClust=3, clMethods=c("pam"), 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' clen_res1 <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
#' 
#' heatmapCmap(clen_res1, "pam", "3", orderMethod="2")
#' heatmapCmap(clen_res1, "pam", "3", orderMethod="2", MultiInstance="concdrug")
#' }
#' @export
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @docType methods
#' @rdname heatmapCmap
#' 
setGeneric("heatmapCmap", 
           function(object, method=clusterMethods(object), 
                    nCluster=nClusters(object), orderMethod="max", MultiInstance="drug",
                    CutoffNumGeneset=20, CutoffPVal=0.05, mergeMethod="mean",
                    low="grey", high="red", na.value="white", maintitle=NULL,
                    printGS=FALSE, add2=TRUE, geom="tile") 
               standardGeneric("heatmapCmap"))

#' @rdname heatmapCmap
#' @aliases heatmapCmap,cogena
setMethod("heatmapCmap", signature(object="cogena"),
          function(object, method=clusterMethods(object), 
                   nCluster=nClusters(object), orderMethod="max", MultiInstance="drug",
                   CutoffNumGeneset=20, CutoffPVal=0.05, mergeMethod="mean",
                   low="grey", high="red", na.value="white", maintitle="cogena",
                   printGS=FALSE, add2=TRUE, geom="tile") {
              
              MultiInstance <- match.arg(MultiInstance, c("drug", "celldrug", "conccelldrug", "concdrug"))
              geom <- match.arg(geom, c("tile", "circle"))
              
              enrichment_score <- enrichment(object, method, nCluster, add2=add2, 
                                             orderMethod=orderMethod, 
                                             CutoffNumGeneset=Inf, 
                                             CutoffPVal=CutoffPVal)
              
              enrichment_score <- as.data.frame(t(enrichment_score))
              
              if (MultiInstance == "drug") {
                  enrichment_score$name <- sapply(strsplit(rownames(enrichment_score) , "@"), "[[", 1)
              } else if (MultiInstance == "celldrug") {
                  enrichment_score$name <- sapply(strsplit(rownames(enrichment_score) , "#"), "[[", 1)
              } else if (MultiInstance == "conccelldrug") {
                  enrichment_score$name <- sapply(strsplit(rownames(enrichment_score) , "_"), "[[", 1)
              } else if (MultiInstance == "concdrug") {
                  enrichment_score$name <- paste(sapply(strsplit(rownames(enrichment_score) , "@|#|_"), "[[", 1), 
                                                 sapply(strsplit(rownames(enrichment_score) , "@|#|_"), "[[", 3), 
                                                 sep="#")
              }
              
              #Multi instance merging
              meanX <- function(x) {
                  instanceCount <- length(which(x>=round(-log2(CutoffPVal))))
                  if (mergeMethod == "max") {
                      suppressWarnings(meanScore <- max(x[which(x>=round(-log2(CutoffPVal)))], na.rm=TRUE) )
                  } else if (mergeMethod == "mean") {
                      meanScore <- mean(x[which(x>=round(-log2(CutoffPVal)))])
                  } else {
                      stop("Parameter mergeMethod should be 'max' or 'min'.")
                  }
                  if (is.nan(meanScore) || is.infinite(meanScore) || instanceCount == 1) {
                      return (0)
                  } else {
                      return (round(meanScore, 1)) 
                  }
              }
              name <- NULL
              score <- dplyr::group_by(enrichment_score, name)
              # score <- dplyr::summarise_each(score, funs(meanX))
              score <- dplyr::summarise_all(score, funs(meanX))
              score <- as.data.frame(score[which(rowSums(score[,-1])!=0), ])
              rownames(score) <- score$name
              if (nrow(score)==0){
                  stop(paste("No enrichment for this order Method:", orderMethod))
              }
              
              score <- t(subset(as.data.frame(score), select=-name))
              score[score==0] <- NA
              
              # the orderMethod options
              # orderMethod <- match.arg(orderMethod, c(rownames(score), "max", "mean"))
              if (orderMethod == "mean") {
                  score = score[,order(colMeans(score, na.rm=TRUE), decreasing=TRUE), drop=FALSE]
                  index_above_cutoffPVal <- which(suppressWarnings(
                      apply(score, 2, max, na.rm=TRUE)) > -log2(CutoffPVal))
              } else if (orderMethod == "max") {
                  colMax <- function(X) {suppressWarnings(
                      apply(X, 2, max, na.rm=TRUE))}
                  score = score[,order(colMax(score), decreasing=TRUE), drop=FALSE]
                  index_above_cutoffPVal <- which(suppressWarnings(
                      apply(score, 2, max, na.rm=TRUE)) > -log2(CutoffPVal))
              } else if (any(grepl( paste0("^", orderMethod, "#"), rownames(score) )) ) {
                  orderMethod_with_num <- rownames(score)[grepl( paste0("^", orderMethod, "#"), rownames(score) )]
                  if (length(orderMethod_with_num) != 1) {
                      stop("Wrong orderMethod!!")
                  }
                  score = score[, order(score[orderMethod_with_num,], decreasing=TRUE), drop=FALSE]
                  index_above_cutoffPVal <- 
                      which(score[orderMethod_with_num,] > -log2(CutoffPVal))
              } else {
                  stop("Wrong orderMethod!")
              }
              
              if (length(index_above_cutoffPVal) > CutoffNumGeneset){
                  score <- score[,c(1:CutoffNumGeneset)]
              } else if (length(index_above_cutoffPVal) == 0){
                  score <- NA
                  return (score)
              } else {
                  #drop para used as length(index_above_cutoffPVal)==1.
                  score <- score[,index_above_cutoffPVal, drop=FALSE]
              }
              
              # Show max CutoffNumGeneset gene sets
              if (nrow(score) > CutoffNumGeneset) {
                  score <- score[,1:CutoffNumGeneset]
              }
              score <- score[,ncol(score):1, drop=FALSE]
              if (printGS==TRUE) {
                  cat (rev(colnames(score)), sep ="\t")
              }
              
              cl_color0 <- upDownGene(object, method, nCluster, add2)
              cl_color <- sapply(cl_color0, function(x) {if (x>=0.6) "red" else if (x<=-0.6) "green4" else "black"})
              cl_color <- c(cl_color, "blue")

              enrich_score <- reshape2::melt(score)
              
              #legend breaks
              cutoff_score <- round(-log2(CutoffPVal), 2)
              max_score <- max(enrich_score$value, na.rm=TRUE)
              if (is.infinite(max_score)) {max_score <- 1000}
              if (max_score/cutoff_score > 2) {
                  breaks <- c(cutoff_score, round( seq(cutoff_score, max_score, length.out=5)[-1]))
              } else if (max_score/cutoff_score >1) {
                  breaks <- c(cutoff_score, round( seq(cutoff_score, max_score, length.out=2)[-1]))
              } else {
                  breaks <- NULL
              }
              
              if (!is.null(maintitle)) {
                  maintitle=paste(maintitle, "\n", "cogena:", method, nCluster)
              } else {
                  maintitle=paste(method, nCluster, "clusters",sep="_")
              }
              
              Var1=Var2=value=NULL

              # ggplot(enrich_score, aes(as.factor(Var1), Var2)) +
              #     geom_tile(aes(fill = value)) +
              #     scale_fill_gradient2("score", space="Lab", mid=low, midpoint=4, low=low,
              #                          high=high, na.value=na.value, breaks=breaks) +
              #     geom_text(aes(fill=value, label=value),size=4, na.rm=TRUE) +
              #     labs(list(title = maintitle, x = "Cluster", y = "Gene set")) +
              #     theme(axis.text.y = element_text(size = rel(1.5), face="bold")) +
              #     theme(axis.text.x = element_text(size = rel(1.3), angle=-90,
              #                                      face="bold", color=cl_color, vjust=0.5))  +
              #     theme(panel.grid.major.x = element_line(color = "grey", size = 5),
              #           panel.grid.major.y = element_blank())
              
              p <- ggplot2::ggplot(enrich_score, aes(as.factor(Var1), Var2)) +
                  labs(list(title = maintitle, x = "Cluster", y = "Gene set")) +
                  theme(axis.text.y = element_text(size = rel(1.5), face="bold")) +
                  theme(axis.text.x = element_text(size = rel(1.3), angle=-90, 
                                                   face="bold", color=cl_color, vjust=0.5))
              
              
              if (geom =="tile") {
                  p +  geom_tile(aes(fill = value)) + 
                      scale_fill_gradient2("score", space="Lab", mid=low, midpoint=4, low=low, 
                                           high=high, na.value=na.value, breaks=breaks) +
                      geom_text(aes(label=value),size=4, na.rm=TRUE)  +
                      theme(panel.grid.major.x = element_line(color = "grey", size = 5),
                            panel.grid.major.y = element_blank())
              } else if (geom =="circle") {
                  p + geom_point(aes(color=value, size = value), na.rm = TRUE, stroke = 3) + 
                      guides(size=FALSE) + 
                      scale_color_gradient2("score", space="Lab", mid=low, midpoint=4, low=low, 
                                            high=high, na.value=na.value, breaks=breaks) +
                      theme(panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey90"),
                            panel.background = element_blank(),
                            axis.line.x = element_line(size = 0.25, linetype = 'solid'),
                            axis.line.y = element_line(size = 0.25, linetype = 'solid'),
                            plot.background = element_blank() )
              }
          }
)

