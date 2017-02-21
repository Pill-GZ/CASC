library(CASC)

load(file = "simul_nonassortative.Rdata")
test_that("casc outputs expected clustering results", {
  case3.result <- CASpecClust(A = adjMat, X = covMat, K = 3, method = 3)
  expect_equal(length(case3.result$cluster), 1500)
  expect_equal(length(case3.result$alpha), 100)
  expect_equal(length(case3.result$WCSS), 100)
})
rm(adjMat, covMat)
