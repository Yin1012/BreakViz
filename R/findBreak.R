#'Filter reads with size requirement
#'
#'@param bedFile A bed file
#'@param minOverlap A length that a part of one read aligned to chromosomes
#'@return the list of read parts met our size requirment.
#'@export
#'@examples
#'bedFile1 <- rtracklayer::import(system.file('extdata', 'test_file_1.bed', package = 'BreakViz'))
#'sizefilter(bedFile1,10)
sizefilter <- function(bedFile, minOverlap = 1000) {
  # filter reads with minimum overlapping size
  SizeFilteredReads <-
    bedFile[width(bedFile@ranges[, ncol(bedFile)]) > minOverlap]
  #throw error if there is no reads met requirment
  if (NROW(SizeFilteredReads) == 0) {
    stop("There is no reads met your minOVerlap requirment")
  }
  return(SizeFilteredReads)
}
#'Select reads that only mapping to two parts of chromosome
#'
#'@param bedFile A bed file after sizeFilter
#'@return the list of read name that only mapping to two parts of chromosome
#'@export
#'@examples
#'bedFile1 <- rtracklayer::import(system.file('extdata', 'test_file_1.bed', package = 'BreakViz'))
#'dupFileter(bedFile1)
dupFileter <- function(bedFile) {
  #search duplication of reads name in all reads
  numOccur <- data.frame(table(bedFile$name))
  preFiltedRead <- numOccur[numOccur$Freq > 1, ]
  #throw error if there is no reads met requirment
  if (NROW(preFiltedRead) == 0) {
    stop("There is no reads mapping into two parts in chromosomes")
  }
  return(preFiltedRead)
}
#'Combine reads that only mapping to two parts of chromosome and two parts are unique to each other
#'
#'@param bedFile A bed file after sizeFilter
#'@param preFiltedRead A file after sizeFilter and dupFilter
#'@return the list of read that only mapping to two parts of chromosome and their information of each part is included
#'@export
#'@examples
#'bedFile2 <- rtracklayer::import(system.file('extdata', 'test_file_2.bed', package = 'BreakViz'), format = 'bed')
#'dupFiltedBedFile2 <- readRDS(system.file('extdata', 'dupFiltedBedFile2.Rda', package = 'BreakViz'))
#'combineSameRead(bedFile2, dupFiltedBedFile2)
combineSameRead <- function(bedFile, preFiltedRead) {
  # create a matrix to store read pair
  prePair <- matrix(NA, NROW(preFiltedRead), 8)
  # create col name in a matrix
  colnames(prePair) <-
    c(
      "read_pair_num",
      "read_name",
      "read_part1_seqnames",
      "read_part1_start",
      "read_part1_width",
      "read_part2_seqnames",
      "read_part2_start",
      "read_part2_width"
    )
  # put two part of one read that have same name together
  n <- 1
  # combine for read part based on each read name
  for (item in preFiltedRead$Var1) {
    if (length(unique(bedFile[which(bedFile$name == item)]@seqnames)) != 1) {
      #set pair number and read name
      prePair[n, 1] <- n
      prePair[n, 2] <- item
      #provide chromosome number of read parts based on sorted read part
      sortedChrom <-
        sortChrom(bedFile[which(bedFile$name == item)])
      #extract first chromosome showed in sortedChrom
      chr1 <- sortedChrom[1, 1]
      numChr1 <-
        as.numeric(gsub("\\D", "", sortedChrom[1, 3]))
      #extract second chromosome showed in sortedChrom
      chr2 <- sortedChrom[2, 1]
      numChr2 <-
        as.numeric(gsub("\\D", "", sortedChrom[2, 3]))
      # extract information based on two extracted read part
      prePair[n, 3] <- chr1
      prePair[n, 4] <-
        bedFile[which(bedFile$name == item)][numChr1]@ranges@start
      prePair[n, 5] <-
        bedFile[which(bedFile$name == item)][numChr1]@ranges@width
      prePair[n, 6] <- chr2
      prePair[n, 7] <-
        bedFile[which(bedFile$name == item)][numChr2]@ranges@start
      prePair[n, 8] <-
        bedFile[which(bedFile$name == item)][numChr2]@ranges@width
      n <- n + 1
    }
  }
  prePair <- prePair[-c(n:NROW(preFiltedRead)),]
  #throw error if there is no read part met requirment
  if (NROW(prePair) == 0) {
    stop("There is no read uniquely mapping into two parts in chromosomes")
  }
  return(prePair)
}
#'Check distance between each parts in one read pair if met the distance requirement
#'
#'@param maxDistance The distance requirment
#'@param read1part1Start The start position of part1 of read1
#'@param read1part1SRange The range of part1 of read1
#'@param read1part2Start The start position of part2 of read1
#'@param read1part2SRange The range of part2 of read1
#'@param read2part1Start The start position of part1 of read2
#'@param read2part1SRange The range of part1 of read2
#'@param read2part2Start The start position of part2 of read2
#'@param read2part2SRange The range of part2 of read2
#'@return The logical value showing if there is a read pair met distance requirment
#'@export
#'@examples
#'#'checkDistance(1000,1,10001,10005,300,2000,200,2230,500)
checkDistanceSingleReadPair <-
  function(maxDistance,
           read1part1Start,
           read1part1Range,
           read1part2Start,
           read1part2Range,
           read2part1Start,
           read2part1Range,
           read2part2Start,
           read2part2Range) {
    # Conver to numeric
    read1part1Start <- as.numeric(gsub("\\D", "", read1part1Start))
    read1part1Range <- as.numeric(gsub("\\D", "", read1part1Range))
    read1part2Start <- as.numeric(gsub("\\D", "", read1part2Start))
    read1part2Range <- as.numeric(gsub("\\D", "", read1part2Range))
    read2part1Start <- as.numeric(gsub("\\D", "", read2part1Start))
    read2part1Range <- as.numeric(gsub("\\D", "", read2part1Range))
    read2part2Start <- as.numeric(gsub("\\D", "", read2part2Start))
    read2part2Range <- as.numeric(gsub("\\D", "", read2part2Range))
    # Calculate the end of read part
    read1part1End <- read1part1Start + read1part1Range
    read1part2End <- read1part2Start + read1part2Range
    read2part1End <- read2part1Start + read2part1Range
    read2part2End <- read2part2Start + read2part2Range
    # there are two conditions. One of part of two read will be first one and distance between two read parts should meet the requirment
    if (((
      read1part1Start - read2part1End < maxDistance &&
      read1part1Start - read2part1End > 0
    ) || (
      read2part1Start -
      read1part1End < maxDistance &&
      read2part1Start - read1part1End > 0
    )
    ) && ((
      read1part2Start - read2part2End <
      maxDistance &&
      read1part2Start - read2part2End > 0
    ) || (
      read2part2Start - read1part2End < maxDistance &&
      read2part2Start - read1part2End > 0
    )
    )) {
      return(TRUE)
    }
    return(FALSE)
  }

#'Check distance requirment on all read pair appeared in prePair list
#'
#'@param prePair read pairs met before quirments
#'@param maxDistance The distance requirment
#'@return A list with possible read pairs
#'@export
#'@examples
#'dupFiltedBedFile2 <- readRDS(system.file('extdata', 'dupFiltedBedFile2.Rda', package = 'BreakViz'))
#'prePair  <- combineSameRead(bedFile2, dupFiltedBedFile2)
#'checkDistanceAllReadPair(prePair, maxDistance)
checkDistanceAllReadPair <- function(prePair, maxDistance) {
  #get number of row in prePair
  numRow <- NROW(prePair)
  # select real pairs
  breakPair <- matrix(NA, numRow, 9)
  colnames(breakPair) <-
    c(
      "read_pair_num",
      "read_name",
      "read_part1_seqnames",
      "read_part1_start",
      "read_part1_width",
      "read_part2_seqnames",
      "read_part2_start",
      "read_part2_width",
      "GeneAnnotation"
    )
  n <- 1
  for (i in 1:numRow) {
    # check if there is breakpoint pair have same chromosome position
    for (item in which(prePair[, 3] == prePair[i, 3])) {
      if (item != i && prePair[item, 6] == prePair[i, 6]) {
        colPair <- 2 * n - 1
        # check the distance between two parts of each read
        if (checkDistanceSingleReadPair(
          maxDistance,
          prePair[item, 4],
          prePair[item, 5],
          prePair[item,
                  7],
          prePair[item, 8],
          prePair[i, 4],
          prePair[i, 5],
          prePair[i, 7],
          prePair[i, 8]
        ) == TRUE) {
          #extract information of read part in prePair
          for (m in 2:NCOL(prePair)) {
            breakPair[colPair, m] <- prePair[item, m]
          }
          for (m in 2:NCOL(prePair)) {
            breakPair[colPair + 1, m] <- prePair[i, m]
          }
          breakPair[colPair, 1] <- n
          breakPair[colPair + 1, 1] <- n
          n <- n + 1
        }
      }
    }
  }
  #only leave useful information
  breakPair <- breakPair[-c(n:100),]
  #throw error if there is mo read pair meeting requirment
  if (NROW(breakPair) == 0) {
    stop("There is no read pair met the min Distance requirment")
  }
  return(breakPair)
}
#'Sort reads parts according to their chromosome number
#'
#'@param lists A line including read part information of one read with chromosome number
#'@return A line whose chromosome number has been ordered.
#'@export
#'@examples
#'bedFile2 <- rtracklayer::import(system.file('extdata', 'test_file_2.bed', package = 'BreakViz'), format = 'bed')
#'dupFiltedBedFile2 <- readRDS(system.file('extdata', 'dupFiltedBedFile2.Rda', package = 'BreakViz'))
#'lists <- bedFile2[which(bedFile2$name == dupFiltedBedFile2$Var1[1])]
#'sortChrom(lists)
sortChrom <- function(lists) {
  #create a new matrix to save sorted chromosomes
  listChrom <- matrix(NA, NROW(lists), 3)
  n <- 1
  # covert chrx and chr Y into number to compare
  for (i in 1:NROW(lists)) {
    item <- lists[i,]
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
  #sort chromosome based on order
  sortedChrom <- as.matrix(listChrom[order(listChrom[, 1]),])
  return(sortedChrom)
}
