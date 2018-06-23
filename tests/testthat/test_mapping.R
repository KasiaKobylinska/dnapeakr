context("mapping")
library(dnapeakr)

test_that("mapping gives empty result for empty input", {
  expect_equal(dnapeakr::mapping(c()), c())
})
