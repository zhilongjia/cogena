% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coExp.R
\name{coExp}
\alias{coExp}
\title{co-expressed gene-set enrichment analysis}
\usage{
coExp(
  obj,
  nClust,
  clMethods = "hierarchical",
  metric = "correlation",
  method = "complete",
  ncore = 2,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{obj}{Differentially expressed gene (DEG) expression profilings. Either 
a numeric matrix, a data.frame, or an ExpressionSet object. Data frames must
contain all numeric columns. In all cases, the rows are the items to be 
clustered (e.g., genes), and the columns are the samples.}

\item{nClust}{A numeric vector giving the numbers of clusters to be 
evaluated. e.g., 2:6 would evaluate the number of clusters ranging from 
2 to 6.}

\item{clMethods}{A character vector giving the clustering methods. The 
default is "hierarchical". Available options are "hierarchical", "kmeans", 
"diana", "fanny", "som", "model", "sota", "pam", "clara", "apcluster", 
and "agnes", with multiple choices allowed.}

\item{metric}{the distance measure to be used. This should be one of 
"euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", 
"abspearson", "correlation", "abscorrelation", "NMI", "biwt", "spearman" or
"kendall". Any unambiguous substring can be given. In detail, please 
reference the parameter method in amap::Dist. Some of the cluster methods
could use only part of the metric. See Detail.}

\item{method}{For hierarchical clustering (hierarchical and agnes), the 
agglomeration method used. The default is "complete". Available choices are
"ward", "single", "complete", and "average".}

\item{ncore}{Number of core used. The default is 2.}

\item{verbose}{verbose.}

\item{...}{to interal function vClusters.}
}
\value{
a genecl object
}
\description{
Co-expressed gene-set enrichment analysis. Gene sets could be 
Pathway, Gene ontology. The gene co-expression is obtained by various
clustering methods.
}
\details{
For metric parameter, "hierarchical","kmeans","diana","fanny","pam" and 
"agnes" can use all the metrics.
"clara" uses "manhattan" or "euclidean", other metric will be changed as 
"euclidean".
"sota" uses "correlation" or "euclidean", other metric will be changed as 
"euclidean".
"model" uses its own metric and "som" and "ap" uses euclidean only, which is 
irrelative with metric.

method:
Available distance measures are (written for two vectors x and y):
\itemize{
\item euclidean Usual square distance between the two vectors (2 norm).
\item maximum Maximum distance between two components of x and y 
(supremum norm).
\item manhattan Absolute distance between the two vectors (1 norm).
\item canberra \eqn{sum(|x_i - y_i| / |x_i + y_i|)} Terms with zero 
numerator and denominator are omitted from the sum and treated as if the 
values were missing.
\item binary (aka asymmetric binary): The vectors are regarded as binary 
bits, so non-zero elements are 'on' and zero elements are 'off'. The 
distance is the proportion of bits in which only one is on amongst those in 
which at least one is on.
\item pearson Also named "not centered Pearson" 
\eqn{1 - sum(x_i y_i) / sqrt [sum(x_i^2) sum(y_i^2)]}.
\item abspearson Absolute Pearson 
\eqn{1 - |sum(x_i y_i) / sqrt [sum(x_i^2) sum(y_i^2)] |}.
\item correlation Also named "Centered Pearson" \eqn{1 - corr(x,y)}.
\item abscorrelation Absolute correlation \eqn{1 - | corr(x,y) |}.
\item spearman Compute a distance based on rank.
\item kendall Compute a distance based on rank. 
\eqn{\sum_{i,j} K_{i,j}(x,y)} with \eqn{K_{i,j}(x,y)} is 0 if \eqn{x_i}, 
\eqn{x_j} in same order as \eqn{y_i}, \eqn{y_j}, 1 if not.
\item NMI normalised mutual information. (use correlation instead so far!)
\item biwt  a weighted correlation based on Tukey's biweight
}
}
\examples{
data(Psoriasis)

#cogena parameters
# the number of clusters. A vector.
nClust <- 2:6
# the number of cores. 
ncore <- 2
# the clustering methods
clMethods <- c("hierarchical","kmeans")
# the distance metric
metric <- "correlation"
# the agglomeration method used for hierarchical clustering (hierarchical 
#and agnes)
method <- "complete"

# See examples of clEnrich function
# the co-expression analysis
\dontrun{
genecl_result <- coExp(DEexprs, nClust=nClust, clMethods=clMethods, 
metric=metric, method=method, ncore=ncore, verbose=TRUE)
}

}
\seealso{
\code{\link{clEnrich}}
}
