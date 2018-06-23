context("appending overlapping genes")
library(dnapeakr)

test_that("overlapping genome data filled for proper peaks", {
  path_znf_714 <- file.path("..", "..", "data", "ZNF714_peaks_processed_score_signal_exo.bed")
  znf_peaks_first_5 <- dnapeakr::mapping(read.table(path_znf_714)[1:5, ], "ZNF714")
  path_small_hg19 <- file.path("..", "..", "data", "hg19_chr1_first_10million.out")
  genome <- dnapeakr::read_human_genomic_file(path_small_hg19)
  appended <- dnapeakr::append_genome_overlaps(znf_peaks_first_5, genome)

  expect_equal(is.na(appended$Start_seq[1]), TRUE)
  expect_equal(is.na(appended$End_seq[1]), TRUE)
  expect_equal(is.na(appended$matching_repeat[1]), TRUE)
  expect_equal(is.na(appended$repeat_class_family[1]), TRUE)
  expect_equal(is.na(appended$length_of_the_overlap[1]), TRUE)

  expect_equal(appended$Start_seq[2], 1440346)
  expect_equal(appended$End_seq[2], 1441337)
  expect_equal(appended$matching_repeat[2], "LTR13")
  expect_equal(appended$repeat_class_family[2], "LTR/ERVK")
  expect_equal(appended$length_of_the_overlap[2], 272)
})
