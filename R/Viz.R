#'Visualize possible pairs
#'
#'@param bedfile A bed file that contains possible pairs for visualization.
#'@param minOverlap  A length that a part of one read aligned to chromosomes
#'@param maxDistance A distance that the parts of one pair between each other
#'@import chromoMap
#'@export
#'@return A list with possible read pairs
#'@examples
#'bedFile <- import(system.file('extdata', 'test_file_4.bed', package = 'BreakViz'), format = 'bed')
#'visPossiblePair(bedFile, minOverlap = 100, maxDistance = 1000, 1)
visPossiblePair <-
  function(bedfile,
           minOverlap = 100,
           maxDistance = 1000,
           colorSet = 1) {
    # convert file format
    newBedfile <- changeFormatChromoMap(bedfile)
    # get four lists including reads after four filters
    sizeFiltedList <- sizefilter(bedfile, minOverlap)
    dupFiltedList <- dupFilter(sizeFiltedList)
    uniqfiltedList <- combineSameRead(bedfile, dupFiltedList)
    distanceFiltedList <-
      checkDistanceAllReadPair(uniqfiltedList, maxDistance)
    # add heat number
    newBedfile <- addHeat(newBedfile, sizeFiltedList$name, 3)
    newBedfile <- addHeat(newBedfile, dupFiltedList$Var1, 3)
    newBedfile <- addHeat(newBedfile, uniqfiltedList[, 2], 5)
    newBedfile <- addHeat(newBedfile, distanceFiltedList[, 2], 15)
    # provide options for color set
    if (colorSet == 1) {
      # green/blue color set
      chromoMap::chromoMap(
        newBedfile,
        type = "heatmap-single",
        chCol = "#EDEDED",
        chBorder = "#EDEDED",
        HeatColRange = c("#A2D145", "#02C39A", "#05668D")
      )
    } else if (colorSet == 2) {
      # red/yellow color set
      chromoMap::chromoMap(
        newBedfile,
        type = "heatmap-single",
        chCol = "#EDEDED",
        chBorder = "#EDEDED",
        HeatColRange = c("#FFD700", "#DC143C")
      )
    } else {
      # yellow/blue color set
      chromoMap::chromoMap(
        newBedfile,
        type = "heatmap-single",
        chCol = "#EDEDED",
        chBorder = "#EDEDED",
        HeatColRange = c("#FFD700", "#1E90FF")
      )
    }
    # provide list to be downloaded if needed
  }
