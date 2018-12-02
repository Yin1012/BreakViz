# BreakViz

## Overview
The goal of BreakViz is to detect possible breakpoints' location in human whole genome sequencing data from MinION and visualize them with by using heatmap.

### Process
After sequencing and mapping, human DNA will be transformed into a bed file.The bed format file will be used as import in my package. User must provides two filters about selecting breakpoints : Mindistance, Maxdistance. We will select possible breakpoints based on these two parameters and create a heat table for whole reads. The package will add heat to reads based on their status after filtering. After visulization, a heatmap version of chromosomes will be showed.

## Dependencies
chromoMap, shiny

## Installation

You can install BreakViz from github with:


``` r
# install.packages("devtools")
devtools::install_github("Yin1012/BreakViz")
```

## Example pipeline

This is a basic example which shows you how to solve a common problem:

``` r


library(devtools)
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
library(rtracklayer)
devtools::load_all(".")

bedFile <- import(system.file("extdata", "test_file_4.bed", package = "BreakViz"), format = "bed")
searchPossiblePairs(bedFile,100, 1000)
visPossiblePair(bedFile, minOverlap = 100, maxDistance = 1000, baseCol = 1)
```
## Sample output

![alt text](http://steipe.biochemistry.utoronto.ca/abc/students/images/b/b8/Result.PNG)
