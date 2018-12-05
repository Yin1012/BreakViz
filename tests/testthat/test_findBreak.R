context("findBreak")
#main functions
# ==== BEGIN SETUP AND PREPARE =================================================
library(rtracklayer)

bedFile1 <- rtracklayer::import(system.file("extdata", "test_file_1.bed", package = "BreakViz"))
bedFile2 <- rtracklayer::import(system.file("extdata", "test_file_2.bed", package = "BreakViz"), format = "bed")
bedFile3 <- rtracklayer::import(system.file("extdata", "test_file_3.bed", package = "BreakViz"), format = "bed")
sizeFiltedBedFile2 <- readRDS(system.file("extdata", "sizeFiltedBedFile2.Rda", package = "BreakViz"))
dupFiltedBedFile2 <- readRDS(system.file("extdata", "dupFiltedBedFile2.Rda", package = "BreakViz"))
dupFiltedBedFile3 <- readRDS(system.file("extdata", "dupFiltedBedFile3.Rda", package = "BreakViz"))
uniqFiltedBedFile2 <- readRDS(system.file("extdata", "uniqFiltedBedFile2.Rda", package = "BreakViz"))

# ==== END SETUP AND PREPARE ===================================================
test_that("There is no read met the size requirment", {
  expect_error(sizefilter(bedFile1,100000), "There is no reads met your minOVerlap requirment")
  expect_error(sizefilter(bedFile2,100000), "There is no reads met your minOVerlap requirment")
  })

test_that("There are reads meet the size requirment", {
  expect_equal(sizefilter(bedFile1,10), bedFile1)
  expect_equal(sizefilter(bedFile2,1000), sizeFiltedBedFile2)
  })

test_that("There is no read showing in two parts in chromosomes", {
  expect_error(dupFileter(bedFile1), "There is no reads mapping into two parts in chromosomes")
                 })

test_that("There are reads showing in two parts in chromosomes", {
  expect_equal(dupFileter(bedFile2), dupFiltedBedFile2)
 })

test_that("There is no read showing in two parts in uniqie chromosomes", {
  expect_error(combineSameRead(bedFile3, dupFiltedBedFile3) ,"There is no read uniquely mapping into two parts in chromosomes")
})

test_that("There are reads showing in two parts in uniqie chromosomes",{
  expect_equal(combineSameRead(bedFile2, dupFiltedBedFile2), uniqFiltedBedFile2)
})

# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]


