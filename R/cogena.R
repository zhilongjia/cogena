#' Co-expressed gene set enrichment analysis
#' 
#' To discovery smaller scale, but highly correlated cellular events 
#' that may be of great biological relevance, co-expressed gene set 
#' enrichment analysis, cogena, clusters gene expression profiles (coExp) and
#' then make enrichment analysis for each clusters (clEnrich) 
#' based on hyper-geometric test. The heatmapCluster and heatmapPEI can
#' visualise the results. See vignette for the detailed workflow.
#'
#' @name cogena_package
#' @aliases cogena cogean_package
#' @docType package
#' @keywords package
#' @source https://github.com/zhilongjia/cogena
#' @examples 
#' 
#' ## A quick start
#' 
#' # Loading the examplar dataseat
#' data(PD)
#' 
#' # Clustering the gene expression profiling
#' clMethods <- c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes")
#' genecl_result <- coExp(DEexprs, nClust=2:6, clMethods=clMethods, 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' # Gene set used
#' annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt", package="cogena")
#' 
#' # Enrichment analysis for clusters
#' clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)
#' 
#' summary(clen_res)
#' 
#' 
#' # Visualisation
#' heatmapCluster(clen_res, "hierarchical", "6")
#' heatmapPEI(clen_res, "hierarchical", "6", printGS=FALSE)
#' 
#' # Obtain genes in a certain cluster
#' head(geneInCluster(clen_res, "hierarchical", "6", "2"))
#' 
#' ## The end
#' 
NULL

#' Parkinson's Disease dataset.
#' 
#' an example dataset of Parkinson's Disease. This dataset is used for 
#' illustration of the usage of \code{cogena} package. It has been normalised
#' the expression profling using \code{rma} method, filtered some 
#' non-informative genes using \code{MetaDE} package and analysed the 
#' differentially expressed genes using \code{limma} package with 
#' the p-value 0.05.
#' 
#' @format three objects: DEexprs, sampleLabel and cogena_result.
#' \describe{
#' \item{DEexprs}{expression of DEG. There are 1243 DEGs and 17 samples.}
#' \item{sampleLabel}{the label of sample, There are 9 control and 8 PD.}
#' \item{cogena_result}{an example of cogena result.}
#' }
#' @docType data
#' @keywords datasets
#' @name PD
#' @aliases PD,DEexprs,sampleLabel,cogena_result
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20163}
NULL

#' gene expression of DEG
#' 
#' @format matrix with 1243 DEGs (row) and 17 samples (column).
#' @docType data
#' @keywords datasets
#' @name DEexprs
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20163}
NULL

#' label of samples
#' 
#' @format a vector with 17 element.
#' @docType data
#' @keywords datasets
#' @name sampleLabel
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20163}
NULL

#' All the gene symbols
#' 
#' @format a vector with 18986 gene symbols.
#' @docType data
#' @keywords datasets
#' @name AllGeneSymbols
#' @source \url{http://www.genenames.org/}
NULL
