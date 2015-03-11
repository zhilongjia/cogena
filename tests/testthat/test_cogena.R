library(cogena)
data(PD)
annoGMT <- "c2.cp.kegg.v4.0.symbols.gmt"
annofile <- system.file("extdata", annoGMT, package="cogena")
anno <- gene2set(annofile, rownames(DEexprs))
data(AllGeneSymbols)
annotationGenesPop <- gene2set(annofile, AllGeneSymbols)
annotationGenesPop <- annotationGenesPop[,colnames(anno)]
nClust <- 2:3
ncore <- 2 #2
clMethods <- c("hierarchical","kmeans")
metric <- "correlation"
method <- "complete"
context("cogena results")

test_that("cogena: metric", {
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric=metric, method=method, annotation=anno, 
        sampleLabel=sampleLabel, ncore=ncore, 
        annotationGenesPop=annotationGenesPop, verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric="euclidean", method=method, annotation=anno, 
        sampleLabel=sampleLabel, ncore=ncore, 
        annotationGenesPop=annotationGenesPop, verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric="spearman", method=method, annotation=anno, 
        sampleLabel=sampleLabel, ncore=ncore, 
        annotationGenesPop=annotationGenesPop, verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric="pearson", method=method, annotation=anno, 
        sampleLabel=sampleLabel, ncore=ncore, 
        annotationGenesPop=annotationGenesPop, verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric="manhattan", method=method, annotation=anno, 
        sampleLabel=sampleLabel, ncore=ncore, 
        annotationGenesPop=annotationGenesPop, verbose=TRUE), "cogena")

})

test_that("cogena: method", {
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric=metric, method=method, annotation=anno, 
        sampleLabel=sampleLabel, ncore=ncore, 
        annotationGenesPop=annotationGenesPop, verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric=metric, method="single", annotation=anno, 
        sampleLabel=sampleLabel, ncore=ncore, 
        annotationGenesPop=annotationGenesPop, verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric=metric, method="average", annotation=anno, 
        sampleLabel=sampleLabel, ncore=ncore, 
        annotationGenesPop=annotationGenesPop, verbose=TRUE), "cogena")
    expect_is(cogena_result <- cogena(DEexprs, nClust=nClust, 
        clMethods=clMethods, metric=metric, method="ward", annotation=anno, 
        sampleLabel=sampleLabel, ncore=ncore, 
        annotationGenesPop=annotationGenesPop, verbose=TRUE), "cogena")

})
