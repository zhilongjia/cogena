library(cogena)
data(PD)
clMethods <- c("hierarchical","kmeans")
genecl_result <- coExp(DEexprs, nClust=2:3, 
                       clMethods=clMethods, 
                       metric="correlation", 
                       method="complete", 
                       ncore=2, 
                       verbose=TRUE)
annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt", package="cogena")
clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)

context("test Class_methods")

test_that("test Class_methods", {
    expect_is(clusterMethods(genecl_result), "character")
    expect_is(clusterMethods(clen_res), "character")

    expect_is(nClusters(genecl_result), "numeric")
    expect_is(nClusters(clen_res), "numeric")

    expect_is(geneclusters(genecl_result, "hierarchical", "3"), "integer")
    expect_is(geneclusters(clen_res, "kmeans", "3"), "integer")

    expect_is(mat(genecl_result), "matrix")
    expect_is(mat(clen_res), "matrix")
    
     
})
