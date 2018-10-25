standChr <- c(1:21, "X", "Y")


#'Visualize possible pairs
#'
#'@param bedfile A bed file
#'@param minOverlap A sizeOverlapping data
#'@param maxDistance The distance requirment
#'@param baseCol The base color
#'@import chromoMap
#'@return A list with possible read pairs
visPossiblePair <- function(bedfile, minOverlap = 1000, maxDistance = 100, baseCol = -10) {
    # convert file format
    newBedfile <- changeFormatChromoMap(bedfile, baseCol)
    # get four lists including reads after four filters
    sizeFiltedList <- sizefilter(bedfile, minOverlap)
    if (identical(sizeFiltedList, "There is no reads met your minOVerlap requirment") == TRUE) {
        return(sizeFiltedList)
    }
    dupFiltedList <- dupFileter(sizeFiltedList)
    if (identical(dupFiltedList, "There is no reads mapping into two parts in chromosomes") == TRUE) {
        return(dupFiltedList)
    }
    uniqfiltedList <- combineSameRead(bedfile, dupFiltedList)
    if (identical(uniqfiltedList, "There is no read uniquely mapping into two parts in chromosomes") ==
        TRUE) {
        return(uniqfiltedList)
    }
    distanceFiltedList <- searchPossiblePairs(bedfile, minOverlap, maxDistance)
    if (identical(distanceFiltedList, "There is no read pair met the min Distance requirment") == TRUE) {
        return("There is no read pair met the min Distance requirment")
    }
    # add heat number
    newBedfile <- addHeat(newBedfile, sizeFiltedList$name, 3)
    newBedfile <- addHeat(newBedfile, dupFiltedList$Var1, 3)
    newBedfile <- addHeat(newBedfile, uniqfiltedList[, 2], 5)
    newBedfile <- addHeat(newBedfile, distanceFiltedList[, 2], 15)
    chromoMap::chromoMap(newBedfile, type = "heatmap-single", chCol = "#EDEDED", chBorder = "#EDEDED", HeatColRange = c("#A2D145", "#02C39A","#05668D"))
    # return (newBedfile)
}



