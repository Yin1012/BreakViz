#'Change format for visualization by chromoMap
#'
#'@param bedfile A bed file
#'@return A bed file filled with base col in each line
#'@export
#'@examples
#'bedFile <- import(system.file('extdata', 'test_file_4.bed', package = 'BreakViz'), format = 'bed')
#'changeFormatChromoMap(bedFile)
changeFormatChromoMap <- function(bedfile) {
    #create a new bedfile
    newBedfile <- as.data.frame(matrix(NA, NROW(bedfile), 4))
    #set col name
    colnames(newBedfile) <- c("name", "chrom", "start", "data")
    #change type to number and integer
    newBedfile$start <- as.numeric(newBedfile$start)
    newBedfile$data <- as.integer(newBedfile$data)
    #set each col based on ChromoMap format
    for (i in 1:NROW(bedfile)) {
        newBedfile[i, 1] <- toString(bedfile[i, ]$name)
        newBedfile[i, 2] <- toString(bedfile[i, ]@seqnames)
        newBedfile[i, 3] <- bedfile[i, ]@ranges@start
        newBedfile[i, 4] <- -10
    }
    return(newBedfile)
}

#'Add heat to the line we want
#'
#'@param formattedBedfile A bed file
#'@param list A list included name of reads we want to add heat
#'@param heatnumber number you want to added
#'@return A bed file after adding heat
#'@export
#'@examples
#'bedFile <- import(system.file('extdata', 'test_file_4.bed', package = 'BreakViz'), format = 'bed')
#'formattedBedfile <- changeFormatChromoMap(bedFile)
#'addHeat(formattedBedfile,1)
addHeat <- function(formattedBedfile, list, heatnumber) {
    for (i in 1:NROW(list)) {
        #select line only appeared in list from formattedBedfile
        for (num in which(formattedBedfile[, 1] == list[i])) {
            formattedBedfile[num, 4] <- as.numeric(formattedBedfile[num, 4]) + as.numeric(heatnumber)
        }
    }
    return(formattedBedfile)
}
