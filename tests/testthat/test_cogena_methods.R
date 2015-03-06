library(cogena)
data(PD)
context("methods of cogena object")

test_that("methods of cogena object", {
    expect_is(clusterMethods(cogena_result), "character")
    expect_is(nClusters(cogena_result), "numeric")
    expect_is(clusters(cogena_result, "kmeans"), "list")
    expect_is(mat(cogena_result), "matrix")
    expect_is(mat(cogena_result), "matrix")
    expect_is(mat(cogena_result), "matrix")
})
