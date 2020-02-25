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
#' data(Psoriasis)
#' 
#' # Clustering the gene expression profiling
#' clMethods <- c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes")
#' genecl_result <- coExp(DEexprs, nClust=5:6, clMethods=clMethods, 
#'     metric="correlation", method="complete", ncore=2, verbose=TRUE)
#' 
#' # Gene set used
#' annofile <- system.file("extdata", "c2.cp.kegg.v7.01.symbols.gmt.xz", package="cogena")
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

#' Psoriasis dataset.
#' 
#' An example dataset of Psoriasis. This dataset is used for 
#' illustration of the usage of \code{cogena} package. It has been normalised
#' the expression profling using \code{rma} method, filtered some 
#' non-informative genes using \code{MetaDE} package and analysed the 
#' differentially expressed genes using \code{limma} package with 
#' the cut-off adjuested p-value 0.05 and abs(logFC) >=1.
#' 
#' @format two objects: DEexprs and sampleLabel.
#' \describe{
#' \item{DEexprs}{expression of DEG. There are 706 DEGs and 116 samples.}
#' \item{sampleLabel}{the label of sample, There are 58 control and 58 Psoriasis.}
#' }
#' @docType data
#' @keywords datasets
#' @name Psoriasis
#' @aliases Psoriasis,DEexprs,sampleLabel,cogena_result
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13355}
NULL

#' gene expression of DEG
#' 
#' @format matrix with 706 DEGs (row) and 116 samples (column).
#' @docType data
#' @keywords datasets
#' @name DEexprs
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13355}
NULL

#' label of samples
#' 
#' @format a vector with 116 element.
#' @docType data
#' @keywords datasets
#' @name sampleLabel
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13355}
NULL

#' All the gene symbols
#' 
#' @format a vector with 18986 gene symbols.
#' @docType data
#' @keywords datasets
#' @name AllGeneSymbols
#' @source \url{http://www.genenames.org/}
NULL
