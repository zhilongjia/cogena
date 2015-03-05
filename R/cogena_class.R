#' An S4 class to represent co-expressed gene-set enrichment analysis result.
#'
#' @slot mat Differentially expressed gene expression profilings. Either a 
#' numeric matrix, a data.frame, or an ExpressionSet object. Data frames must
#' contain all numeric columns. In all cases, the rows are the items to be 
#' clustered (e.g., genes), and the columns are the samples.
#' @slot clusterObjs a list contains clustering results.
#' @slot Distmat the distance matrix.
#' @slot measures a list of the enrichment results.
#' @slot clMethods clustering method.
#' @slot labels the label of genes
#' @slot nClust A numeric vector giving the numbers of clusters to be 
#' evaluated. e.g., 2:6 would evaluate the number of clusters ranging from 2 to 6.
#' @slot metric the distance measure to be used. This must be one of "euclidean",
#' "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", 
#' "correlation", "abscorrelation", "spearman" or "kendall". Any unambiguous 
#' substring can be given. In detail, please reference the parameter method in 
#' amap::Dist. Some of the cluster methods could use only part of the metric. 
#' Please reference the manual of cogena. 
#' @slot method For hierarchical clustering (hclust and agnes), the agglomeration 
#' method used. The default is "complete". Available choices are "ward", "single", 
#' "complete", and "average".
#' @slot annotation logical matrix of biological annotation with row be DE gene, 
#' column be gene sets and value be logical. 
#' @slot sampleLabel character vector with names are sample names. only used for plotting.
#' @slot ncore the number of cores used.
#' @slot call the called function
#' @rdname cogena_class
#' @exportClass cogena
#' @import class


#setClassUnion("character or factor", c("character", "factor"))
setClass("cogena", slots= list(mat="matrix",
                          clusterObjs="list",
                          Distmat="dist",
                          measures="list",
                          clMethods="character",
                          labels="character",
                          nClust="numeric",
                          metric="character", 
                          method="character",
                          annotation="matrix", 
                          sampleLabel="factor",
                          ncore="numeric",
                          call="call"))

