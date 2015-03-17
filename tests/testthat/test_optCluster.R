library(cogena)
data(PD)
context("test optCluster function")

test_that("test based parameter in optCluster function", {
    expect_is(optCluster(cogena_result), "matrix")
    expect_is(optCluster(cogena_result, based="I"), "matrix")
})
