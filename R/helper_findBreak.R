#'Sort reads parts according to their chromosome number
#'
#'@param lists A line including read part information of one read with chromosome number
#'@return A line whose chromosome number has been ordered.
sortChrom <- function(lists) {
  listChrom <- matrix(NA, NROW(lists), 3)
  n <- 1
  # covert chrx and chr Y into number to compare
  for (i in 1:NROW(lists)) {
    item <- lists[i, ]
    if (toString(item@seqnames) == "chrX") {
      listChrom[n, 1] <- 22
      listChrom[n, 2] <- "chrX"
    } else if (toString(item@seqnames) == "chrY") {
      listChrom[n, 1] <- 23
      listChrom[n, 2] <- "chrY"
    } else {
      listChrom[n, 1] <- as.numeric(gsub("\\D", "", item@seqnames))
      listChrom[n, 2] <- toString(item@seqnames)
    }
    listChrom[n, 3] <- n
    n <- n + 1
  }
  sortedChrom <- as.matrix(listChrom[order(listChrom[, 1]), ])
  return(sortedChrom)
}

#'Copy information from outFile to inFile
#'
#'@param inFile File we want to transfer into
#'@param numRowInFile The number of rows in inFile we want to transfer information
#'@param outFile File we want to transfer from
#'@param numRowOutFile The number of rows in inFile we want to transfer information into
copyInfo <- function(inFile, numRowInFile, outFile, numRowOutFile, numCol) {

  for (i in 2:numCol) {

    inFile[numRowInFile, i] <- outFile[numRowOutFile, i]
  }
}
