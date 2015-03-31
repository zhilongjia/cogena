#' co-expressed gene-set enrichment analysis
#' 
#' Co-expressed gene-set enrichment analysis. Gene sets could be 
#' Pathway, Gene ontology. The gene co-expression is obtained by various
#' clustering methods.
#' 
#' For metric parameter, "hierarchical","kmeans","diana","fanny","pam" and 
#' "agnes" can use all the metrics.
#' "clara" uses "manhattan" or "euclidean", other metric will be changed as 
#' "euclidean".
#' "sota" uses "correlation" or "euclidean", other metric will be changed as 
#' "euclidean".
#' "model" uses its own metric and "som" uses euclidean only, which is 
#' irrelative with metric.
#'
#' method:
#' Available distance measures are (written for two vectors x and y):
#' \itemize{
#' \item euclidean Usual square distance between the two vectors (2 norm).
#' \item maximum Maximum distance between two components of x and y 
#' (supremum norm).
#' \item manhattan Absolute distance between the two vectors (1 norm).
#' \item canberra \eqn{sum(|x_i - y_i| / |x_i + y_i|)} Terms with zero 
#' numerator and denominator are omitted from the sum and treated as if the 
#' values were missing.
#' \item binary (aka asymmetric binary): The vectors are regarded as binary 
#' bits, so non-zero elements are 'on' and zero elements are 'off'. The 
#' distance is the proportion of bits in which only one is on amongst those in
#'  which at least one is on.
#' \item pearson Also named "not centered Pearson" 
#' \eqn{1 - sum(x_i y_i) / sqrt [sum(x_i^2) sum(y_i^2)]}.
#' \item abspearson Absolute Pearson 
#' \eqn{1 - |sum(x_i y_i) / sqrt [sum(x_i^2) sum(y_i^2)] |}.
#' \item correlation Also named "Centered Pearson" \eqn{1 - corr(x,y)}.
#' \item abscorrelation Absolute correlation \eqn{1 - | corr(x,y) |}.
#' \item spearman Compute a distance based on rank.
#' \item kendall Compute a distance based on rank. 
#' \eqn{\sum_{i,j} K_{i,j}(x,y)} with \eqn{K_{i,j}(x,y)} is 0 if \eqn{x_i}, 
#' \eqn{x_j} in same order as \eqn{y_i}, \eqn{y_j}, 1 if not.
#' \item NMI normalised mutual information. (use correlation instead so far!)
#' \item biwt  a weighted correlation based on Tukey's biweight
#' }
#' 
#' @param obj Differentially expressed gene (DEG) expression profilings. Either 
#' a numeric matrix, a data.frame, or an ExpressionSet object. Data frames must
#' contain all numeric columns. In all cases, the rows are the items to be 
#' clustered (e.g., genes), and the columns are the samples.
#' @param nClust A numeric vector giving the numbers of clusters to be 
#' evaluated. e.g., 2:6 would evaluate the number of clusters ranging from 
#' 2 to 6.
#' @param clMethods A character vector giving the clustering methods. The 
#' default is "hierarchical". Available options are "hierarchical", "kmeans", 
#' "diana", "fanny", "som", "model", "sota", "pam", "clara", and "agnes", 
#' with multiple choices allowed.
#' @param metric the distance measure to be used. This should be one of 
#' "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", 
#' "abspearson", "correlation", "abscorrelation", "NMI", "biwt", "spearman" or
#'  "kendall". Any unambiguous 
#' substring can be given. In detail, please reference the parameter method in 
#' amap::Dist. Some of the cluster methods could use only part of the metric. 
#' See Detail.
#' @param method For hierarchical clustering (hierarchical and agnes), the 
#' agglomeration method used. The default is "complete". Available choices are
#'  "ward", "single", "complete", and "average".
#' @param annofile gene set filename.
#' @param sampleLabel factor or character vector with names are sample names. 
#' only used for plotting.
#' @param ncore Number of core used. The default is 2.
#' @param TermFreq a value from [0,1) to filter low-frequence gene sets.
#' @param verbose verbose.
#' @param ... to interal function vClusters.
#' @return a cogena object
#' @examples
#' data(PD)
#' 
#' #annotaion
#' annoGMT <- "c2.cp.kegg.v4.0.symbols.gmt"
#' annofile <- system.file("extdata", annoGMT, package="cogena")
#' #cogena parameters
#' # the number of clusters. A vector.
#' nClust <- 2:6
#' # the number of cores. 
#' ncore <- 2
#' # the clustering methods
#' clMethods <- c("hierarchical","kmeans")
#' # the distance metric
#' metric <- "correlation"
#' # the agglomeration method used for hierarchical clustering (hierarchical 
#' #and agnes)
#' method <- "complete"
#' 
#' # the cogena analysis
#' cogena_result <- cogena(DEexprs, nClust=nClust, clMethods=clMethods, 
#'     metric=metric, method=method,  annofile=annofile, sampleLabel=sampleLabel, 
#'     ncore=ncore, verbose=TRUE)
#'
#' @export
#' @import devtools

cogena <- function(obj, nClust, clMethods="hierarchical",
                    metric="correlation", method="complete", 
                    annofile=NULL, sampleLabel=NULL, ncore=2, 
                    TermFreq=0, verbose=FALSE,...) {

  clMethods <- match.arg(clMethods, c("hierarchical","kmeans","diana","fanny",
                                     "som","model","sota","pam","clara",
                                     "agnes"), several.ok=TRUE)
  
  ## used for hierarchical, kmeans, diana, fanny, agnes, pam
  metric <- match.arg(metric,c("euclidean", "correlation", "abscorrelation",
   "manhattan", "spearman", "maximum", "kendall", "canberra", "binary", 
   "pearson", "abspearson", "NMI", "biwt")) #
  
#   if("NMI" %in% metric){
# #       if(!require(infotheo)) {
# #           stop("package 'infotheo' required for mutual information metric")
#       devtools::install_github("zhilongjia/infotheo")
#       require(infotheo)
#   }
  
  
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
         stop("argument 'obj' must be a matrix, data.frame, or ExpressionSet 
              object")
         )
  
  if(is.null(rownames(mat))) {
      stop("rownames of data must be present")
  }
  
  if (is.null(annofile)) {
      annofile <- system.file("extdata", "c2.cp.kegg.v4.0.symbols.gmt", package="cogena")
  }
  
  annotation <- gene2set(annofile, rownames(mat), TermFreq=TermFreq)
  # the background gene gene-sets matrix
  annotationGenesPop <- gene2set(annofile, AllGeneSymbols, TermFreq=TermFreq)
  annotationGenesPop <- annotationGenesPop[,colnames(annotation)]

  if (ncol(annotationGenesPop) ==0 || ncol(annotation)==0) {
    stop("Error in annotation or annotationGenesPop as ncol equals zero. 
         Maybe lower the TermFreq value.")
  }

  

  nClust <- floor(nClust)
  if (any(nClust<1))
      {stop("argument 'nClust' must be a positive integer vector")}

  #Dist caculation
  if (verbose) {print("Dist caculation")}
  if (metric == "biwt") {
        Distmat <- as.dist(biwt::biwt.cor(mat, output="distance"))
    } else if (metric == "NMI") {
        #NMI from infotheo package
        #nmi <- infotheo::nMI(infotheo::discretize(t(mat)), method= "emp")
        #Distmat <- as.dist(1-nmi)
        Distmat <- amap::Dist(mat, method="correlation", nbproc=ncore)
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
    if (verbose) print(paste("##### The clMethod,", clMethods[i], "starts ##"))
    
    cvalid <- vClusters(mat, Distmat , clMethods[i], nClust, method=method, 
                        metric=metric, annotation=annotation,
                        ncore=ncore, annotationGenesPop, verbose=verbose, ...)
    clusterObjs[[i]] <- cvalid$clusterObj
    validMeasures[[i]] <- cvalid$measures
    if (verbose) print(paste("The clMethod,", clMethods[i], "done"))
  }
  

  if (is.character(sampleLabel)){
      sampleLabel <- as.factor(sampleLabel)
  }
  res <- new("cogena", mat=mat, Distmat=Distmat, clusterObjs=clusterObjs, 
      measures=validMeasures, clMethods=clMethods, labels=rownames(mat), 
      nClust=nClust, metric=metric, method=method, annotation=annotation, 
      sampleLabel=sampleLabel, ncore=ncore, gmt=basename(annofile), call=match.call())
  return (res)
}
