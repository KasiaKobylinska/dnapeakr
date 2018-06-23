context("mapping")
library(dnapeakr)

test_that("mapping gives empty result for empty input", {
  expect_equal(dnapeakr::mapping(c()), c())
})

test_that("some values are found for first row of ZNF714", {
  path_znf_714 <- file.path("..", "..", "data", "ZNF714_peaks_processed_score_signal_exo.bed")
  first_row_znf714 <- read.table(path_znf_714)[1, ]
  result_first_row <- dnapeakr::mapping(first_row_znf714)
  expect_equal(is.na(result_first_row$ensembl_gene_id), FALSE)
  expect_equal(result_first_row$chromosome_name, 1)
  expect_equal(is.na(result_first_row$start_position), FALSE)
  expect_equal(is.na(result_first_row$end_position), FALSE)
})
