library(cogena)
data(PD)
context("test enrichment function")

test_that("test orderMethod parameter in enrichment function", {
    expect_is(enrichment(cogena_result, "kmeans", "3", orderMethod = "I"), "matrix")
    expect_is(enrichment(cogena_result, "kmeans", "3", orderMethod = "2"), "matrix")
    expect_is(enrichment(cogena_result, "kmeans", "3", orderMethod = "II"), "logical")
})
