# reduce bed file for MCscan to longest transcript

#####!/usr/bin/env Rscript


library(tidyverse)
infile <- "data/Anc0.bed"

#args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
#  stop("At least one argument must be supplied (input file).bed", call.=FALSE)
#} else if (length(args)==1) {
#  # default output file
#  infile = args[0]
#}

#infile <- "data/Anc0.bed"

#setwd("/media/kate/Massey_linux_onl/projects/Ancestral/MCscan/")
bed <- read.delim(infile, header = FALSE, stringsAsFactors = FALSE)
bed$length <- abs(bed$V3-bed$V2)
bed <- bed %>% group_by(V4) %>% top_n(1, length)
bed$length <- NULL

write.table(bed, sub(".bed", ".maxlen.bed", infile), row.names = FALSE, quote = FALSE)


# reformat and subset to regions for samtools faidx
bed <- bed %>% tidyr::unite("pos", V2:V3, sep = "-", remove = FALSE)
bed <- bed %>% tidyr::unite("faidx", V1:pos, sep = ":", remove = FALSE)
#bed <- bed %>% mutate(pos = paste(V2, V3, collapse = "-"),
#                      faidx = paste(V1, pos, collapse = ":"))

bed <- bed[, c("faidx", "V4")]
#bed <- as.matrix(bed)
#colnames(bed) <- NULL
head(bed)

write.table(bed, sub(".bed", "_regions_lookup.txt", infile), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

bed <- bed[, c("faidx")]
write.table(bed, sub(".bed", "_regions.txt", infile), row.names = FALSE, quote = FALSE, col.names = FALSE)


