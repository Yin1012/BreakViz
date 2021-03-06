% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_viz.R
\name{addHeat}
\alias{addHeat}
\title{Add heat to the line we want}
\usage{
addHeat(formattedBedfile, list, heatnumber)
}
\arguments{
\item{formattedBedfile}{A bed file}

\item{list}{A list included name of reads we want to add heat}

\item{heatnumber}{number you want to added}
}
\value{
A bed file after adding heat
}
\description{
Add heat to the line we want
}
\examples{
bedFile <- import(system.file('extdata', 'test_file_4.bed', package = 'BreakViz'), format = 'bed')
formattedBedfile <- changeFormatChromoMap(bedFile)
addHeat(formattedBedfile,1)
}
