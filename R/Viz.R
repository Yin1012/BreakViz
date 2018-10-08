standChr <- c(1:21, "X","Y")

changeFormatChromoMap <- function(bedfile, baseCol){
  newBedfile <-  as.data.frame(matrix(NA, NROW(bedfile), 4))
  colnames(newBedfile) <- c("name", "chrom" , "start", "data")
  newBedfile$start <- as.numeric(newBedfile$start)
  newBedfile$data <- as.integer(newBedfile$data)
  for (i in 1:NROW(bedfile)){
    newBedfile[i, 1] <- toString(bedfile[i,]$name)
    newBedfile[i, 2] <- toString(bedfile[i,]@seqnames)
    newBedfile[i, 3] <- bedfile[i,]@ranges@start
    newBedfile[i, 4] <- baseCol
  }
  return (newBedfile)
  }

addHeat <- function(formattedBedfile, list, heatnumber){
  for (i in 1: NROW(list)){
    for (num in which(formattedBedfile[,1] == list[i])){
      formattedBedfile[num, 4] <- as.numeric(formattedBedfile[num, 4]) + as.numeric(heatnumber)
    }
  }

  return (formattedBedfile)
}

visPossiblePair <- function(bedfile,minOverlap = 1000, maxDistance = 100, baseCol = 1){
  #convert file format
  newBedfile <- changeFormatChromoMap(bedfile, baseCol)
  #get four lists including reads after four filters
  sizeFiltedList <- sizefilter(bedfile, minOverlap)
  if (identical(sizeFiltedList, "There is no reads met your minOVerlap requirment") == TRUE){
    return (sizeFiltedList)
  }
  dupFiltedList <- dupFileter(sizeFiltedList)
  if (identical(dupFiltedList,"There is no reads mapping into two parts in chromosomes") == TRUE){
    return (dupFiltedList)
  }
  uniqfiltedList <- combineSameRead(bedfile, dupFiltedList)
  if (identical(uniqfiltedList, "There is no read uniquely mapping into two parts in chromosomes") == TRUE) {
    return (uniqfiltedList)
  }
  distanceFiltedList <- searchPossiblePairs(bedfile, minOverlap, maxDistance)
  if (identical(distanceFiltedList,"There is no read pair met the min Distance requirment") == TRUE) {
    return ("There is no read pair met the min Distance requirment")
  }
  #add heat number
  newBedfile <- addHeat(newBedfile, sizeFiltedList$name, 1)
  newBedfile <- addHeat(newBedfile, dupFiltedList$Var1, 1)
  newBedfile <- addHeat(newBedfile, uniqfiltedList[,2], 1)
  newBedfile <- addHeat(newBedfile, distanceFiltedList[,2], 1)
  chromoMap::chromoMap(newBedfile, type = "heatmap-single")
  return (newBedfile)
}

visReads <- function(bedfile){
  newBedfile <- changeFormatChromoMap(bedfile, baseCol)
  chromoMap::chromoMap(newBedfile, type = "annotation")
}

showBreakPair <- function(prePairs, prePairsNum, read1Col, read2Col, figCols){
  prePairs <- as.data.frame(prePairs)
  pair <- as.data.frame(matrix(NA, 2, 5))
  colnames(pair) <- c("Chrom", "Start", "End", "Name", "Colors")
  pair$Start <- as.integer(pair$Start)
  pair$End <- as.integer(pair$End)
  pair$Chrom <- as.character(pair$Chrom)
  chr <- c()
  n <- 1
  for (item in which(prePairs[,1] == prePairsNum)){
    pair$Chrom[n] <- paste("chr", prePairs[item, 3], sep = "")
    pair$Name[n] <- prePair[item, 2]
    pair$Start[n] <- as.integer(toString(prePairs[item, 4]))
    print(pair$Start[n])
    pair$End[n] <- as.integer(toString(prePairs[item, 4])) + as.integer(toString(prePairs[item,5]))
    pair$Colors[n] <- read1Col
    pair$Chrom[n + 1] <- paste("chr", prePairs[item, 6], sep = "")
    pair$Name[n + 1] <- prePair[item, 2]
    pair$Start[n + 1] <- as.integer(toString(prePairs[item, 7]))
    print(pair$Start[n+1])
    pair$End[n + 1] <- as.integer(toString(prePairs[item, 7])) + as.integer(toString(prePairs[item,8]))
    pair$Colors[n + 1] <- read2Col
    chr[1] <- prePairs[item, 3]
    chr[2] <- prePairs[item, 6]
    n <- n + 1
  }
  return (pair)
  chromPlot::chromPlot(gaps=hg_gap, bands=pair, chr, figCols)
  return (pair)
}


