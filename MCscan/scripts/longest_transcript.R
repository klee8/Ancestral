# reduce bed file for MCscan to longest transcript

#####!/usr/bin/env Rscript


library(tidyverse)

#infile <- "data/cactus_Anc0.bed"

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).bed", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  infile = args[1] #| "data/cactus_Anc0.bed"
  infile
}

# read in the bed file, calculate length of transcripts
bed <- read.delim(infile, header = FALSE, stringsAsFactors = FALSE)
bed$length <- abs(bed$V3-bed$V2)
bed <- bed %>% group_by(V4) %>% top_n(1, length)
bed$length <- NULL

write.table(bed, sub(".bed", ".maxlen.bed", infile), row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)


# reformat and subset to regions for samtools faidx
bed <- bed %>% tidyr::unite("pos", V2:V3, sep = "-", remove = FALSE)
bed <- bed %>% tidyr::unite("faidx", V1:pos, sep = ":", remove = FALSE)

# create lookup table
bed <- bed[, c("V4", "faidx")]
write.table(bed, sub(".bed", ".regions_lookup.txt", infile), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

# create file of regions for samtools faidx
bed <- bed[, c("faidx")]
write.table(bed, sub(".bed", ".regions.txt", infile), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")


