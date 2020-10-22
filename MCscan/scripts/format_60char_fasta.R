# change fasta format to 60 chars wide for CHRONicle

library(seqinr)

# get the fasta file name from the command line
args <- commandArgs(trailingOnly = TRUE)
print(args)
input <- args[1]

# read in as a named list
fasta <- read.fasta(input, forceDNAtolower = FALSE)

# add '*' to the end of each protein seq
add_star <- function(seq_vec) {
  if (tail(vector, n=1) != "*") { seq_vec <- append(seq_vec, "*") }
  return(seq_vec)
}
#fasta <- sapply(fasta, function(x) add_star(x))

# write out in new format
write.fasta(fasta, names(fasta), paste("fmt", input, sep = "_") , open = "w", nbchar = 60, as.string = FALSE)
