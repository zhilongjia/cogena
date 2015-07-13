#' An S4 class to represent co-expressed gene
#'
#' @slot mat Differentially expressed gene expression profilings. Either a 
#' numeric matrix, a data.frame, or an ExpressionSet object. Data frames must
#' contain all numeric columns. In all cases, the rows are the items to be 
#' clustered (e.g., genes), and the columns are the samples.
#' @slot clusterObjs a list contains clustering results.
#' @slot Distmat the distance matrix.
#' @slot clMethods clustering method.
#' @slot labels the label of genes
#' @slot nClust A numeric vector giving the numbers of clusters to be 
#' evaluated. e.g., 2:6 would evaluate the number of clusters ranging from 
#' 2 to 6.
#' @slot metric the distance measure to be used. It must be one of 
#' "euclidean","maximum", "manhattan", "canberra", "binary", 
#' "pearson", "abspearson", "correlation", "abscorrelation", 
#' "spearman" or "kendall". Any unambiguous substring can be 
#' given. In detail, please reference the parameter method in 
#' amap::Dist. Some of the cluster methods could use only part 
#' of the metric. Please reference the manual of cogena. 
#' @slot method For hierarchical clustering (hclust and agnes), the 
#' agglomeration method used. The default is "complete". Available 
#' choices are "ward", "single", "complete", and "average".
#' @slot ncore the number of cores used.
#' @slot call the called function
#' @rdname genecl_class
#' @exportClass genecl
#' @import class

setClass("genecl", slots=list(mat="matrix", 
                              clusterObjs="list", 
                              Distmat="dist",
                              clMethods="character",
                              labels="character",
                              nClust="numeric",
                              metric="character", 
                              method="character",
                              ncore="numeric",
                              call="call"))
