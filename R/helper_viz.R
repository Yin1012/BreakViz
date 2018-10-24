#'Change format for visualization by chromoMap
#'
#'@param bedfile A bed file
#'@param baseCol basecolour parameter
#'@return A bed file filled with base col in each line
changeFormatChromoMap <- function(bedfile, baseCol) {
  newBedfile <- as.data.frame(matrix(NA, NROW(bedfile), 4))
  colnames(newBedfile) <- c("name", "chrom", "start", "data")
  newBedfile$start <- as.numeric(newBedfile$start)
  newBedfile$data <- as.integer(newBedfile$data)
  for (i in 1:NROW(bedfile)) {
    newBedfile[i, 1] <- toString(bedfile[i, ]$name)
    newBedfile[i, 2] <- toString(bedfile[i, ]@seqnames)
    newBedfile[i, 3] <- bedfile[i, ]@ranges@start
    newBedfile[i, 4] <- baseCol
  }
  return(newBedfile)
}

#'Add heat to the line we want
#'
#'@param formattedBedfile A bed file
#'@param list A list included name of reads we want to add heat
#'@param heatnumber number you want to added
#'@return A bed file filled with base col in each line
addHeat <- function(formattedBedfile, list, heatnumber) {
  for (i in 1:NROW(list)) {
    for (num in which(formattedBedfile[, 1] == list[i])) {
      formattedBedfile[num, 4] <- as.numeric(formattedBedfile[num, 4]) + as.numeric(heatnumber)
    }
  }

  return(formattedBedfile)
}
