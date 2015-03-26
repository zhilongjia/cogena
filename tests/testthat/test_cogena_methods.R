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

context("methods of cogena object")

test_that("methods of cogena object", {
    expect_is(clusterMethods(cogena_result), "character")
    expect_is(nClusters(cogena_result), "numeric")
    expect_is(clusters(cogena_result, "kmeans"), "list")
})
