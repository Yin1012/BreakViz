#performance.R exdata dev
#


#scenario:
# function to split string into codons
# input e.g "GCATTCTAGGATCAA"
# output: "GCA". "TTC", "TAG", "GAT", "CAA"

S <- "GCATTCTAGGATCAA"
nchar(S)

substring(S,c(1,4), c(3,6))

codon_1 <- function(x) {
  1 <-char(x)
  first <- seq(1,1, by =3)
  last <- first + 2
  return(substring(x,first, last))
}

condon_2 <- function(x){
unilist(strsplit(S, "(>,,,,)", perl = True))
}

codon_3 <- function(x){
  return(unilist(stri_extract_all_regax(x, "...")))
}
#============================
#
x <- paste(sample(c("A", "C", "G", "T"), 99, replace = TRUE), collapse = "")
#============================
Sys.time()
as.integer(Sys.time())

start <- Sys.time()
output <- codon_1(x)
end <- Sys.time()
end - start

system.time({output <- codons_1(x)})

microbenchmark(output <- codons_1(x))

