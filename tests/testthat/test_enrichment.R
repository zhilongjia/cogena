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
context("test enrichment function")

test_that("test orderMethod parameter in enrichment function", {
    expect_is(enrichment(cogena_result, "kmeans", "3", orderMethod = "I"), "matrix")
    expect_is(enrichment(cogena_result, "kmeans", "3", orderMethod = "2"), "matrix")
    expect_is(enrichment(cogena_result, "kmeans", "3", orderMethod = "II"), "logical")
})
