# merge Progressive Cactus and MGLO output using gff gene locations to 
# bridge the two datasets 
# AIM: use gene order from MGLO to set order of progressive cactus blocks



library(tidyverse)



#' Gene order supplied by MGLO in fasta-style format and re-formated to one
#' line prior to reading into R (3 columns: genome, chr, gene_list)
#' Here split the 3rd column (with the list of gene numbers) into a vector
#' and return it named with "genome_chr"
get_gene_order <- function(x) {
  genome = x[1]
  chr = x[2]
  listname <- paste(genome, chr, sep = "_")
  gene_list = x[3]
  gene_order_vec <- gene_list %>% 
    str_trim( side = c("both")) %>% 
    strsplit( ",") %>%
    as.vector()
  names(gene_order_vec) <- listname
  return(gene_order_vec)
}


#####  Get gene locations
# get gene list from the MGLO (genomes at the tips of the tree)
# compare against gff files to get their physical locations

# read in gene order files of tip genomes
tip_gen <- read.delim("data/gene_order_fmtR.dat", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(tip_gen) <- c("genome", "chr", "gene_order")

# create a list of vectors of gene order named by genome and chr ("genome_chr").
tip_gene_lists <- apply(tip_gen, 1, get_gene_order)
str(tip_gene_lists)

# read in gff files


# get gene locations for each genome and merge with 'numeric' gene names from MGLO


# check that list of genes in each block does not clash with order in genome order vectors from MGLO



##### Apply gene locations to maf

# read in re-formatted cactus maf file  (from `perl fmt_cactus_for_R.pl cactus.maf`)
maf <- read.delim("data/cactus.fmtR.maf", stringsAsFactors = FALSE, header = FALSE)
colnames(maf) <- c("block", "s", "genome.chr", "pos", "length", "strand", "chrlength", "sequence")
# remove s column
maf$s <- NULL
# remove blocks with less than 500bp
submaf <- maf %>% group_by(block) %>% filter(!any(mean(length) < 500))

# add list of genes to each block


# check number of blocks without genes





##### get ancestral gene lists and order blocks

# read in gene order files of ancestral genomes
anc_gen <- read.delim("data/geneorder_fmtR.out", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(anc_gen) <- c("genome", "chr", "gene_order")
# create a list of vectors of gene order named by genome and chr ("genome_chr").
anc_gene_lists <- apply(tip_gen, 1, get_gene_order)
str(gene_lists)

# order blocks for each maf 'ancestral' genome by ancestral gene list














