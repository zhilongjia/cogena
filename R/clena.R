#' cluster enrichment analysis
#' 
#' clena : Gene sets enrichment analysis for the cluster of gene expression profilings.
#' Gene sets could be pathway, GO etc.
#' 
#' For metric parameter, "hierarchical","kmeans","diana","fanny","pam" and "agnes" can use
#' all the metrics.
#' "clara” uses "manhattan" or "euclidean", other metric will be changed as "euclidean".
#' "sota" uses "correlation" or "euclidean", other metric will be changed as "euclidean".
#' "model" uses its own metric and "som" uses euclidean only, which is irrelative with metric.
#'
#' method:
#' Available distance measures are (written for two vectors x and y):
#' \itemize{
#' \item euclidean Usual square distance between the two vectors (2 norm).
#' \item maximum Maximum distance between two components of x and y (supremum norm).
#' \item manhattan Absolute distance between the two vectors (1 norm).
#' \item canberra sum(|x_i - y_i| / |x_i + y_i|) Terms with zero numerator and denominator are omitted from the sum and treated as if the values were missing.
#' \item binary (aka asymmetric binary): The vectors are regarded as binary bits, so non-zero elements are ‘on’ and zero elements are ‘off’. The distance is the proportion of bits in which only one is on amongst those in which at least one is on.
#' \item pearson Also named "not centered Pearson" 1 - sum(x_i y_i) / sqrt [sum(x_i^2) sum(y_i^2)].
#' \item abspearson Absolute Pearson 1 - |sum(x_i y_i) / sqrt [sum(x_i^2) sum(y_i^2)] |.
#' \item correlation Also named "Centered Pearson" 1 - corr(x,y).
#' \item abscorrelation Absolute correlation 1 - | corr(x,y) |.
#' \item spearman Compute a distance based on rank.
#' \item kendall Compute a distance based on rank. ∑_{i,j} K_{i,j}(x,y) with K_{i,j}(x,y) is 0 if x_i, x_j in same order as y_i,y_j, 1 if not.
#' }
#' 
#' @param obj Differentially expressed gene expression profilings. Either a 
#' numeric matrix, a data.frame, or an ExpressionSet object. Data frames must
#' contain all numeric columns. In all cases, the rows are the items to be 
#' clustered (e.g., genes), and the columns are the samples.
#' @param nClust A numeric vector giving the numbers of clusters to be 
#' evaluated. e.g., 4:6 would evaluate the number of clusters ranging from 4 to 6.
#' @param clMethods A character vector giving the clustering methods. The 
#' default is "hierarchical". Available options are "hierarchical", "kmeans", 
#' "diana", "fanny", "som", "model", "sota", "pam", "clara", and "agnes", 
#' with multiple choices allowed.
#' @param metric the distance measure to be used. This must be one of "euclidean",
#' "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", 
#' "correlation", "abscorrelation", "spearman" or "kendall". Any unambiguous 
#' substring can be given. In detail, please reference the parameter method in 
#' amap::Dist. Some of the cluster methods could use only part of the metric.
#' @param method For hierarchical clustering (hierarchical and agnes), the agglomeration 
#' method used. The default is "complete". Available choices are "ward", "single", 
#' "complete", and "average".
#' @param annotation logical matrix of biological annotation with row be DE gene, 
#' column be gene sets and value be logical. 
#' @param sampleLabel character vector with names are sample names. only used for plotting.
#' @param ncore Number of core used. The default is 2.
#' @param annotationGenesPop logical matrix of biological annotation with row be 
#' all genes, column be gene sets and value be logical.
#' @param verbose verbose.
#' @return a clena object
#' @examples
#' data("PD")
#' 
#' #annotaion
#' annoGMT <- "c2.cp.v4.0.symbols.gmt"
#' #annoGMT <- "~/clena/data-raw/c2.cp.biocarta.v4.0.symbols.gmt"
#' anno = gene2set(anno=annoGMT, difSAM_GSE7621_DEG, TermFreq=0)
#' annotationGenesPop = gene2set(anno=annoGMT, rownames(GSE7621.filtered.expr), TermFreq=0)
#' annotationGenesPop <- annotationGenesPop[,colnames(anno)]
#' 
#' #clena parameters
#' nClust <- 2:4
#' ncore <- 3
#' clMethods <- c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes")
#'
#' GSE7621.clena.cluster <- clena(GSE7621.SAM.DEG.expr, 
#'                               nClust=nClust, 
#'                               clMethods=clMethods, 
#'                               metric="spearman", 
#'                               method="complete",  
#'                               annotation=anno,  
#'                               sampleLabel=sampleLabel.GSE7621, 
#'                               ncore=ncore, 
#'                               annotationGenesPop=annotationGenesPop, 
#'                               verbose=TRUE)
#'
#' @export

clena <- function(obj, nClust, clMethods="hierarchical",
                    metric="correlation", method="complete", 
                    annotation=NULL, sampleLabel=NULL, ncore=2, 
                    annotationGenesPop=NULL, verbose=FALSE,...) {

  clMethods <- tolower(clMethods)
  clMethods <- match.arg(clMethods,c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes"), several.ok=TRUE)

  if("som"%in%clMethods) {
    if(!require(kohonen)) {
      stop("package 'kohonen' required for clustering using SOM")
    }
  }

  if("model"%in%clMethods) {
    if(!require(mclust)) {
      stop("package 'mclust' required for model-based clustering")
    }
  }
  
  ## used for hierarchical, kmeans, diana, fanny, agnes, pam
  metric <- match.arg(metric,c("euclidean", "correlation", "abscorrelation", "manhattan", 
  "spearman", "maximum", "kendall", "canberra", "binary", "pearson", "abspearson")) 
  
  ## for hclust, agnes
  method <- match.arg(method,c("ward", "single", "complete", "average")) 

  switch(class(obj),
         matrix = mat <- obj,
         ExpressionSet = mat <- Biobase::exprs(obj),
         data.frame = {
           if(any(!sapply(obj,class)%in%c("numeric","integer")))
             stop("data frame 'obj' contains non-numeric data")
           mat <- as.matrix(obj)
         },
         stop("argument 'obj' must be a matrix, data.frame, or ExpressionSet object")
         )

  if ("clara"%in%clMethods & metric %in% c("euclidean", "manhattan"))
    warning("'clara' currently only works with 'euclidean' or 'manhattan' metrics - other metrics will be changed to 'euclidean'  ")

  if (is.null(annotation) || is.null(annotationGenesPop)) {
    stop("annotation must be specified")
  }

  if (ncol(annotationGenesPop) ==0 || ncol(annotation)==0) {
    stop("Error in annotation or annotationGenesPop as ncol equals zero. Maybe lower the TermFreq value.")
  }

  if(is.null(rownames(mat))) {
    stop("rownames of data must be present")
  }

  nClust <- floor(nClust)
  if (any(nClust<1))
    stop("argument 'nClust' must be a positive integer vector")

  ##################################################################
  
  clusterObjs <- vector("list",length(clMethods))
  names(clusterObjs) <- clMethods

  validMeasures = vector("list", length=length(clMethods))
  names(validMeasures) <- clMethods
  
  #for each clMethods, implement the Cluster Analysis.
  for (i in 1:length(clMethods)) {
    if (verbose) print(paste("##### The clMethod,", clMethods[i], "starts #####"))
    
    cvalid <- vClusters(mat, clMethods[i], nClust, method=method, 
                        metric=metric, annotation=annotation,
                      ncore=ncore, annotationGenesPop, verbose=verbose)
    clusterObjs[[i]] <- cvalid$clusterObj
    validMeasures[[i]] <- cvalid$measures
    if (verbose) print(paste("The clMethod,", clMethods[i], "done"))
  }
  
  if(is.null(rownames(mat))) {
    rownames(mat) <- 1:nrow(mat)
    warning("rownames for data not specified, using 1:nrow(data)")
  }

  new("clena", mat=mat, clusterObjs=clusterObjs, measures=validMeasures,
      clMethods=clMethods, labels=rownames(mat), nClust=nClust, 
      metric=metric, method=method, annotation=annotation, 
      sampleLabel=sampleLabel, ncore=ncore, call=match.call())
}


#' clusterMethods : Methods of cluster.
#' @param object a clena object
setGeneric("clusterMethods", function(object, ...) standardGeneric("clusterMethods"))
#' @method clusterMethods
#' @rdname clena
#' @exportMethod clusterMethods
setMethod("clusterMethods",signature(object="clena"),
          function(object) return(object@clMethods))


#' get the number of clusters.
#' @inheritParams clusterMethods
## number of clusters accessor
setGeneric("nClusters", function(object, ...) standardGeneric("nClusters"))
#' @method nClusters
#' @rdname clena
#' @exportMethod nClusters
setMethod("nClusters",signature(object="clena"),
          function(object) return(object@nClust))


#' get the cluster of certain clustering method.
#' @inheritParams clusterMethods
#' @param method as clMethods in clena function
## clusters accessor
setGeneric("clusters", function(object, method, ...) standardGeneric("clusters"))
#' @method clusters
#' @rdname clena
#' @exportMethod clusters
setMethod("clusters",signature(object="clena"),
          function(object, method=clusterMethods(object)) {
              method <- match.arg(method, clusterMethods(object))
              return(object@clusterObjs[[method]])})


#' get the original data from a clena object.
#' @inheritParams clusterMethods
#' @return a matrix
setGeneric("mat", function(object, ...) standardGeneric("mat"))
#' @method : mat
#' @rdname clena
#' @exportMethod mat
setMethod("mat",signature(object="clena"),
          function(object) return(object@mat))


#' get the enrichment table from a clena object.
#' @inheritParams clusters
#' @param nClusters as nClust in clena function.
setGeneric("enrichment", function(object, method, nClusters, ...) standardGeneric("enrichment"))
#' @method enrichment
#' @rdname clena
#' @exportMethod enrichment
setMethod("enrichment",signature(object="clena"),
          function(object,method=clusterMethods(object), nClusters=nClusters(object)) {
              method <- match.arg(method, clusterMethods(object))
              nClusters <- match.arg(nClusters, as.character(nClusters(object)))
              enrichment <- object@measures[[method]][[nClusters]]
              colnames(enrichment) <- tolower(colnames(enrichment))
              enrichment <- enrichment[rownames(enrichment),]
              #enrichment <- enrichment[order(rownames(enrichment)),]
              return (enrichment)
          })

#' summary : a summary of clena object.
#' @method summary
#' @rdname clena
#' @exportMethod summary
setMethod("summary","clena",
          function(object, digits = max(3,getOption("digits")-3)) {
              cat("\nClustering Methods:\n",clusterMethods(object),"\n\n")
              cat("The Number of Clusters:\n",nClusters(object),"\n\n")
              cat("Metric of Distance Matrix:\n", object@metric, "\n\n")
              cat("Agglomeration method for hierarchical clustering (hclust and agnes):\n", object@method, "\n\n")
          })


