# BreakViz

The goal of BreakViz is to detect possible breakpoints' location in human whole genome sequencing data from MinION

## Installation

You can install BreakViz from github with:


``` r
# install.packages("devtools")
devtools::install_github("Yin1012/BreakViz")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```
library(devtools)
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
library(rtracklayer)
devtools::load_all(".")

bedFile <- import(system.file("extdata", "test_file_4.bed", package = "BreakViz"), format = "bed")
searchPossiblePairs(bedFile1,1000, 100)
