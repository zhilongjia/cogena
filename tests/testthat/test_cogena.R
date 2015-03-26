library(cogena)
data(PD)
annoGMT <- "c2.cp.kegg.v4.0.symbols.gmt"
annofile <- system.file("extdata", annoGMT, package="cogena")
nClust <- 2:3
ncore <- 2 #2
clMethods <- c("hierarchical","kmeans")
metric <- "correlation"
method <- "complete"
context("cogena results")

test_that("cogena: metric", {
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric=metric, method=method, annofile=annofile, 
        sampleLabel=sampleLabel, ncore=ncore, 
        verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric="euclidean", method=method, annofile=annofile, 
        sampleLabel=sampleLabel, ncore=ncore, 
        verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric="spearman", method=method, annofile=annofile, 
        sampleLabel=sampleLabel, ncore=ncore, 
        verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric="pearson", method=method, annofile=annofile, 
        sampleLabel=sampleLabel, ncore=ncore, 
        verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric="manhattan", method=method, annofile=annofile, 
        sampleLabel=sampleLabel, ncore=ncore, 
        verbose=TRUE), "cogena")

})

test_that("cogena: method", {
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric=metric, method=method, annofile=annofile, 
        sampleLabel=sampleLabel, ncore=ncore, 
        verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric=metric, method="single", annofile=annofile, 
        sampleLabel=sampleLabel, ncore=ncore, 
        verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric=metric, method="average", annofile=annofile, 
        sampleLabel=sampleLabel, ncore=ncore, 
        verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric=metric, method="ward", annofile=annofile, 
        sampleLabel=sampleLabel, ncore=ncore, 
        verbose=TRUE), "cogena")

})
