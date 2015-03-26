library(cogena)
data(PD)

annoGMT <- "c2.cp.kegg.v4.0.symbols.gmt"
annofile <- system.file("extdata", annoGMT, package="cogena")
nClust <- 2:3
ncore <- 1
clMethods <- c("hierarchical","kmeans")
metric <- "correlation"
method <- "complete"
cogena_result <- cogena(DEexprs, nClust=nClust, clMethods=clMethods, 
    metric=metric, method=method,  annofile=annofile, sampleLabel=sampleLabel, 
    ncore=ncore, verbose=TRUE)

context("test optCluster function")

test_that("test based parameter in optCluster function", {
    expect_is(optCluster(cogena_result, ncores=1), "matrix")
    expect_is(optCluster(cogena_result, based="I", ncores=1), "matrix")
})
