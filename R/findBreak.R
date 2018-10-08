#'Filt reads with size requirment
#'
#'@param bedfile A bed file
#'@param minOverlap A sizeOverlapping data
#'@return the list of read parts met our size requirment.
sizefilter <- function(bedfile,  minOverlap = 1000){
  # filter reads with minimum overlapping size
  SizeFilteredReads <- bedfile[width(bedfile@ranges[,ncol(bedfile)]) > minOverlap]
  if (NROW(SizeFilteredReads) == 0 ){
  return ("There is no reads met your minOVerlap requirment")
  }
  return (SizeFilteredReads)
}

#'Selet reads that only mapping to two parts of chromosome
#'
#'@param bedfile A bed file after sizeFilter
#'@return the list of read name that only mapping to two parts of chromosome
dupFileter <- function(file){
  #search duplication of reads name in all reads
  numOccur <- data.frame(table(file$name))
  preFiltedRead <- numOccur[numOccur$Freq > 1,]
  if (NROW(preFiltedRead) == 0 ){
    return ("There is no reads mapping into two parts in chromosomes")
  }
  return (preFiltedRead)
}

#'Combine reads that only mapping to two parts of chromosome and two parts are unique to each other
#'
#'@param bedfile A bed file after sizeFilter
#'@param preFiltedRead A file after sizeFilter and dupFilter
#'@return the list of read that only mapping to two parts of chromosome and their information of each part is included
combineSameRead <- function(bedfile,preFiltedRead){
  #create a matrix to store read pair
  prePair <- matrix(NA, NROW(preFiltedRead),8)
  colnames(prePair) <- c("read_pair_num", "read_name", "read_part1_seqnames", "read_part1_start", "read_part1_width", "read_part2_seqnames", "read_part2_start", "read_part2_width")
  #put two part of one read that have same name together
  n <- 1
  for (item in preFiltedRead$Var1){
    if (length(unique(bedfile[which(bedfile$name == item)]@seqnames)) != 1) {
    prePair[n, 1] <- n
    prePair[n, 2] <- item
    sortedChrom <- sortChrom(bedfile[which(bedfile$name == item)])
    chr1 <- sortedChrom[1, 1]
    numChr1 <- as.numeric(gsub("\\D","",sortedChrom[1, 3]))
    chr2 <- sortedChrom[2, 1]
    numChr2 <- as.numeric(gsub("\\D","",sortedChrom[2, 3]))
    #print(numChr2)
      prePair[n, 3] <- chr1
      prePair[n, 4] <- bedfile[which(bedfile$name == item)][numChr1]@ranges@start
      prePair[n, 5] <- bedfile[which(bedfile$name == item)][numChr1]@ranges@width
      prePair[n, 6] <- chr2
      prePair[n, 7] <- bedfile[which(bedfile$name == item)][numChr2]@ranges@start
      prePair[n, 8] <- bedfile[which(bedfile$name == item)][numChr2]@ranges@width
    n <- n + 1
    }}
  if (NROW(prePair) == 1 && is.na(prePair[1,1])){
    return ("There is no read uniquely mapping into two parts in chromosomes")
  }
  prePair <- prePair[-c(n:NROW(preFiltedRead)),]
  return (prePair)
}

#'Sort reads parts according to their chromosome number
#'
#'@param lists A line including read part information of one read with chromosome number
#'@return A line whose chromosome number has been ordered.
sortChrom <- function(lists){
  listChrom <- matrix(NA, NROW(lists), 3)
  n <- 1
  #covert chrx and chr Y into number to compare
  for (i in 1: NROW(lists)){
    item <- lists[i,]
    if (toString(item@seqnames) == "chrX"){
      listChrom[n, 1] <- 22
      listChrom[n, 2] <- "chrX"
    } else if (toString(item@seqnames) == "chrY"){
      listChrom[n, 1] <- 23
      listChrom[n, 2] <- "chrY"
    } else {
      listChrom[n, 1] <- as.numeric(gsub("\\D","", item@seqnames))
      listChrom[n, 2] <- toString(item@seqnames)
    }
   listChrom[n, 3] <- n
   n <- n + 1
  }
  sortedChrom <- as.matrix(listChrom[order(listChrom[,1]),])
  return (sortedChrom)
}

#'Copy information from outFile to inFile
#'
#'@param inFile File we want to transfer into
#'@param numRowInFile The number of rows in inFile we want to transfer information
#'@param outFile File we want to transfer from
#'@param numRowOutFile The number of rows in inFile we want to transfer information into
copyInfo <- function(inFile, numRowInFile, outFile, numRowOutFile, numCol){

  for (i in 2:numCol) {

    inFile[numRowInFile, i] <- outFile[numRowOutFile, i]
    }
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
checkDistanceReadPair <- function(maxDistance,read1part1Start, read1part1Range, read1part2Start, read1part2Range, read2part1Start, read2part1Range, read2part2Start, read2part2Range){
  #Conver to numeric
  read1part1Start <- as.numeric(gsub("\\D","",read1part1Start))
  read1part1Range <- as.numeric(gsub("\\D","",read1part1Range))
  read1part2Start <- as.numeric(gsub("\\D","",read1part2Start))
  read1part2Range <- as.numeric(gsub("\\D","",read1part2Range))
  read2part1Start <- as.numeric(gsub("\\D","",read2part1Start))
  read2part1Range <- as.numeric(gsub("\\D","",read2part1Range))
  read2part2Start <- as.numeric(gsub("\\D","",read2part2Start))
  read2part2Range <- as.numeric(gsub("\\D","",read2part2Range))
  #Calculate the end of read part
  read1part1End <- read1part1Start + read1part1Range
  read1part2End <- read1part2Start + read1part2Range
  read2part1End <- read2part1Start + read2part1Range
  read2part2End <- read2part2Start + read2part2Range

  if (((read1part1Start - read2part1End < maxDistance && read1part1Start - read2part1End > 0)|| (read2part1Start - read1part1End < maxDistance && read2part1Start - read1part1End > 0)) &&
    ((read1part2Start - read2part2End < maxDistance && read1part2Start - read2part2End > 0)|| (read2part2Start -read1part2End < maxDistance && read2part2Start -read1part2End > 0))){
    return (TRUE)
  }
  return (FALSE)
}

#'Search possible pairs
#'
#'@param bedfile A bed file
#'@param minOverlap A sizeOverlapping data
#'@param maxDistance The distance requirment
#'@return A list with possible read pairs
searchPossiblePairs <- function(bedfile,  minOverlap = 1000, maxDistance = 100){
  #Filted reads that have occured in two places and meet our overlap requirement
  SizeFilteredReads <- sizefilter(bedfile, minOverlap)
  if (identical(SizeFilteredReads, "There is no reads met your minOVerlap requirment") == TRUE){
    return (SizeFilteredReads)
  }
  preFiltedRead <- dupFileter(SizeFilteredReads)
  if (identical(preFiltedRead,"There is no reads mapping into two parts in chromosomes") == TRUE){
    return ("There is no reads mapping into two parts in chromosomes")
  }
  #Combine reads' informatin together for next step selection
  prePair <- combineSameRead(bedfile, preFiltedRead)
  if (identical(prePair, "There is no read uniquely mapping into two parts in chromosomes") == TRUE) {
    return ("There is no read uniquely mapping into two parts in chromosomes")
  }
  #select real pairs
  breakPair <- matrix(NA, 100, 9)
  colnames(breakPair) <- c("read_pair_num", "read_name", "read_part1_seqnames", "read_part1_start", "read_part1_width", "read_part2_seqnames", "read_part2_start", "read_part2_width", "GeneAnnotation")
  numRow<- NROW(prePair)
  n <- 1
  for (i in 1:numRow){
      #check if there is breakpoint pair have same chromosome position
      for (item in which(prePair[, 3] == prePair[i, 3])) {
      if (item != i && prePair[item, 6] == prePair[i, 6] ){
        colPair <- 2*n - 1
        #check the distance between two parts of each read
        if (checkDistanceReadPair(maxDistance, prePair[item, 4] , prePair[item, 5], prePair[item, 7], prePair[item, 8], prePair[i, 4], prePair[i, 5], prePair[i, 7], prePair[i, 8]) == TRUE){
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
  breakPair <- breakPair[-c(n:100),]
  if (NROW(breakPair) == 0) {
    return ("There is no read pair met the min Distance requirment")
  }
  return (breakPair)
}
  # if two place of reads are found, check if the read size is suitble
  #  to our minOverlap
  # create a table to store data read pair that met requirment
  #

#'Add geneannotation to one read
#'
#'@param readSeqnames An integer representing readSeqnames
#'@param readStart A position showing read start
#'@param readWidth An integer showing read width
#'@param dataType datatype relating to a database
#'@return A list with possible read pairs
addGeneAnnotation <- function(readSeqnames, readStart, readWidth, dataType = "hg38"){
  if (dataType == "hg38"){
    #load("ens.gene.ann.hg38.rda")
    data <- ens.gene.ann.hg38
  } else if (dataType == "hg18"){
    #load("~/BreakViz/data/ens.gene.ann.hg18.rda")
    data <- ens.gene.ann.hg18
  } else if (dataType == "hg19"){
    #load("~/BreakViz/data/ens.gene.ann.hg19.rda")
    data <- ens.gene.ann.hg10
  } else {
    return ("There is no data that meet your requirement")
  }
  readStart <- as.numeric(gsub("\\D","",readStart))
  readWidth <- as.numeric(gsub("\\D","",readWidth))
  geneAnnotationLeft <- data[intersect(intersect(which(data$Start <= readStart), which(data$End >= readStart)), which(data$Chromosome == readSeqnames)),1]
  geneAnnotationRight <- data[intersect(intersect(which(data$Start <= readStart + readWidth), which(data$End >= readStart + readWidth)), which(data$Chromosome == readSeqnames)),1]
  geneAnnotationInside <- data[intersect(intersect(which(data$Start >= readStart), which(data$End <= readStart + readWidth)), which(data$Chromosome == readSeqnames)),1]
  geneAnnotation <- union(union(geneAnnotationLeft,geneAnnotationRight),geneAnnotationInside)
  return (geneAnnotation)
}

#'Add geneannotation to each read in a file
#'
#'@param breakPair An file including read pairs
#'@param dataType datatype relating to a database
#'@return A list with read pairs that have been added gene annotation
getGeneAnnotation <- function(breakPair, dataType = "hg38"){
  for (i in 1: NROW(breakPair)){
    numChr1 <- as.numeric(gsub("\\D","",breakPair[i, 3]))
    numChr2 <- as.numeric(gsub("\\D","",breakPair[i, 6]))
    breakPair[i, 9] <- paste(union(addGeneAnnotation(numChr1, breakPair[i, 4], breakPair[i, 5], dataType), addGeneAnnotation(numChr2, breakPair[i, 7], breakPair[i, 8], dataType)), collapse = " ")
    return (breakPair)
    }
}

