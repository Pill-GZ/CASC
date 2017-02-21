library(CASC)

load("test_casc.Rda")
test_that("casc outputs expected clustering results", {
  case1.result <- CASpecClust(A = A, X = X, K = K)
  expect_equal(length(case1.result$cluster), 10)
  expect_equal(length(case1.result$alpha), 100)
  expect_equal(length(case1.result$WCSS), 100)
})
rm(A, X, K)
