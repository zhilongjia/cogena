library(cogena)
data(PD)
clMethods <- c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes")
clMethods <- c("hierarchical","kmeans","diana","fanny","som","sota","pam","clara","agnes")

annofile <- system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt", package="cogena")

context("Testing coExp, clEnrich and optCluster\n")

test_that("coExp, clEnrich and optCluster", {
	expect_is(genecl_result <- coExp(DEexprs, nClust=2:3, clMethods=clMethods, metric="correlation", method="complete", ncore=2, verbose=TRUE), "genecl")
    expect_is(clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=2),
              "cogena")
    expect_is(optCluster(clen_res, ncores=2), "matrix")

})

