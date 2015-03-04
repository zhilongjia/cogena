#' Parkinson's Disease dataset.
#' 
#' an example dataset of Parkinson's Disease. This dataset is used for 
#' illustration of the usage of \code{cogena} package. It has been normalised
#' the expression profling using \code{rma} method, filtered some non-informative
#' genes using \code{MetaDE} package and analysed the differentially expressed
#' genes using \code{limma} package with the p-value 0.05.
#' 
#' @format three objects: DEexprs, sampleLabel and cogena_result.
#' \describe{
#'  \item{DEexprs}{expression of DEG. There are 1243 DEGs and 17 samples.}
#'  \item{sampleLabel}{the label of sample, There are 9 control and 8 PD.}
#'  \item{cogena_result}{an example of cogena result.}
#' }
#' @docType data
#' @keywords datasets
#' @name PD
#' @aliases PD,DEexprs,sampleLabel,cogena_result
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20163}
NULL
