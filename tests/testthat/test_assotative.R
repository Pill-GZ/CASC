library(CASC)

load(file = "simul_assortative.Rdata")
test_that("casc outputs expected clustering results", {
  case2.result <- CASpecClust(A = adjMat, X = covMat, K = 3, method = 2)
  expect_equal(length(case2.result$cluster), 1500)
  expect_equal(length(case2.result$alpha), 100)
  expect_equal(length(case2.result$WCSS), 100)
})
rm(adjMat, covMat)

