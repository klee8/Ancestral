# read in multiple alignment format file (maf) from cactus exports
setwd("/media/kate/Massey_linux_onl/projects/Ancestral/merge_MGLO_cactus/")

library(tidyverse)

# read in re-formatted cactus maf file  (from `perl fmt_cactus_for_R.pl cactus.maf`)
maf <- read.delim("cactus.fmtR.maf", stringsAsFactors = FALSE, header = FALSE)
colnames(maf) <- c("block", "s", "genome.chr", "pos", "length", "strand", "chrlength", "sequence")
# remove s column
maf$s <- NULL
# remove blocks with less than 500bp
submaf <- maf %>% group_by(block) %>% filter(!any(mean(length) < 500))

# read in gene order files of tip genomes
tipgen <- read.delim("gene_order.dat", header = FALSE, stringsAsFactors = FALSE, sep = ">")
head(tipgen)

