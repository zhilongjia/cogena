#' co-expressed gene-set enrichment analysis
#' 
#' Co-expressed gene-set enrichment analysis. Gene sets could be 
#' Pathway, Gene ontology. The gene co-expression is obtained by various
#' clustering methods.
#' 
#' For metric parameter, "hierarchical","kmeans","diana","fanny","pam" and 
#' "agnes" can use all the metrics.
#' "clara" uses "manhattan" or "euclidean", other metric will be changed as "euclidean".
#' "sota" uses "correlation" or "euclidean", other metric will be changed as "euclidean".
#' "model" uses its own metric and "som" uses euclidean only, which is irrelative with metric.
#'
#' method:
#' Available distance measures are (written for two vectors x and y):
#' \itemize{
#' \item euclidean Usual square distance between the two vectors (2 norm).
#' \item maximum Maximum distance between two components of x and y (supremum norm).
#' \item manhattan Absolute distance between the two vectors (1 norm).
#' \item canberra \eqn{sum(|x_i - y_i| / |x_i + y_i|)} Terms with zero numerator and denominator are omitted from the sum and treated as if the values were missing.
#' \item binary (aka asymmetric binary): The vectors are regarded as binary bits, so non-zero elements are 'on' and zero elements are 'off'. The distance is the proportion of bits in which only one is on amongst those in which at least one is on.
#' \item pearson Also named "not centered Pearson" \eqn{1 - sum(x_i y_i) / sqrt [sum(x_i^2) sum(y_i^2)]}.
#' \item abspearson Absolute Pearson \eqn{1 - |sum(x_i y_i) / sqrt [sum(x_i^2) sum(y_i^2)] |}.
#' \item correlation Also named "Centered Pearson" \eqn{1 - corr(x,y)}.
#' \item abscorrelation Absolute correlation \eqn{1 - | corr(x,y) |}.
#' \item spearman Compute a distance based on rank.
#' \item kendall Compute a distance based on rank. \eqn{\sum_{i,j} K_{i,j}(x,y)} with \eqn{K_{i,j}(x,y)} is 0 if \eqn{x_i}, \eqn{x_j} in same order as \eqn{y_i}, \eqn{y_j}, 1 if not.
#' \item NMI normalised mutual information.
#' \item biwt  a weighted correlation based on Tukey's biweight
#' }
#' 
#' @param obj Differentially expressed gene (DEG) expression profilings. Either a 
#' numeric matrix, a data.frame, or an ExpressionSet object. Data frames must
#' contain all numeric columns. In all cases, the rows are the items to be 
#' clustered (e.g., genes), and the columns are the samples.
#' @param nClust A numeric vector giving the numbers of clusters to be 
#' evaluated. e.g., 4:6 would evaluate the number of clusters ranging from 4 to 6.
#' @param clMethods A character vector giving the clustering methods. The 
#' default is "hierarchical". Available options are "hierarchical", "kmeans", 
#' "diana", "fanny", "som", "model", "sota", "pam", "clara", and "agnes", 
#' with multiple choices allowed.
#' @param metric the distance measure to be used. This should be one of "euclidean",
#' "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", 
#' "correlation", "abscorrelation", "NMI", "biwt", "spearman" or "kendall". Any unambiguous 
#' substring can be given. In detail, please reference the parameter method in 
#' amap::Dist. Some of the cluster methods could use only part of the metric. See Detail.
#' @param method For hierarchical clustering (hierarchical and agnes), the agglomeration 
#' method used. The default is "complete". Available choices are "ward", "single", 
#' "complete", and "average".
#' @param annotation logical matrix of biological annotation with row be DEG, 
#' column be gene sets and value be logical. 
#' @param sampleLabel factor vector with names are sample names. only used for plotting.
#' @param ncore Number of core used. The default is 2.
#' @param annotationGenesPop logical matrix of biological annotation with row be 
#' all genes, column be gene sets and value be logical.
#' @param verbose verbose.
#' @param ... to function vClusters.
#' @return a cogena object
#' @examples
#' data("PD")
#' 
#' #annotaion
#' annoGMT <- "c2.cp.kegg.v4.0.symbols.gmt"
#' annofile <- system.file("data", annoGMT, package="cogena")
#' # the DEG gene-sets matrix
#' anno <- gene2set(annofile, rownames(DEexprs))
#' # the background gene gene-sets matrix
#' annotationGenesPop <- gene2set(annofile, BGgenes)
#' annotationGenesPop <- annotationGenesPop[,colnames(anno)]
#' 
#' #cogena parameters
#' # the number of clusters. A vector.
#' nClust <- 2:3
#' # the number of cores. 
#' ncore <- 2
#' # the clustering methods
#' clMethods <- c("hierarchical","kmeans")
#' # the distance metric
#' metric <- "correlation"
#' # the agglomeration method used for hierarchical clustering (hierarchical and agnes)
#' method <- "complete"
#' 
#' # the cogena analysis
#' cogena_result <- cogena(DEexprs, nClust=nClust, clMethods=clMethods, metric=metric,
#'                       method=method,  annotation=anno, sampleLabel=sampleLabel, 
#'                       ncore=ncore, annotationGenesPop=annotationGenesPop, verbose=TRUE)
#'
#' @export

cogena <- function(obj, nClust, clMethods="hierarchical",
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
  "spearman", "maximum", "kendall", "canberra", "binary", "pearson", "abspearson",
  "NMI", "biwt"))
  
  if("NMI" %in% metric){
#       if(!require(infotheo)) {
#           stop("package 'infotheo' required for mutual information metric")
      devtools::install_github("zhilongjia/infotheo")
      require(infotheo)
  }
  
  if ("biwt" %in% metric){
      if(!require(biwt)) {
          stop("package 'biwt' required for weighted correlation based on Tukey’s biweight")
      }
  }
  
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

  if ("clara" %in% clMethods & !(metric %in% c("euclidean", "manhattan")))
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

    #Dist caculation
    if (verbose) {print("Dist caculation")}
    if (metric == "biwt") {
        Distmat <- as.dist(biwt.cor(mat, output="distance"))
    } else if (metric == "NMI") {
        #NMI from infotheo package
        nmi <- infotheo::NMI(infotheo::discretize(t(mat)), method= "emp")
        Distmat <- as.dist(1-nmi)
        #Distmat <- as.dist(1- (nmi-min(nmi))/(max(nmi)-min(nmi)) )
    } else {
    Distmat <- amap::Dist(mat, method=metric, nbproc=ncore)
    }
    #convert the NA to max of the Distmat
    #Distmat[which(is.na(Distmat))]= max(Distmat, na.rm=T)
    if (verbose) {print("Dist caculation done")}

  ##################################################################
  
  clusterObjs <- vector("list",length(clMethods))
  names(clusterObjs) <- clMethods

  validMeasures = vector("list", length=length(clMethods))
  names(validMeasures) <- clMethods
  
  #for each clMethods, implement the Cluster Analysis.
  for (i in 1:length(clMethods)) {
    if (verbose) print(paste("##### The clMethod,", clMethods[i], "starts #####"))
    
    cvalid <- vClusters(mat, Distmat , clMethods[i], nClust, method=method, 
                        metric=metric, annotation=annotation,
                        ncore=ncore, annotationGenesPop, verbose=verbose, ...)
    clusterObjs[[i]] <- cvalid$clusterObj
    validMeasures[[i]] <- cvalid$measures
    if (verbose) print(paste("The clMethod,", clMethods[i], "done"))
  }
  
  if(is.null(rownames(mat))) {
    rownames(mat) <- 1:nrow(mat)
    warning("rownames for data not specified, using 1:nrow(data)")
  }
  if (is.character(sampleLabel)){
      sampleLabel <- as.factor(sampleLabel)
  }
  new("cogena", mat=mat, Distmat=Distmat, clusterObjs=clusterObjs, measures=validMeasures,
      clMethods=clMethods, labels=rownames(mat), nClust=nClust, 
      metric=metric, method=method, annotation=annotation, 
      sampleLabel=sampleLabel, ncore=ncore, call=match.call())
}




