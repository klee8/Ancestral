# merge Progressive Cactus and MGLO output using gff gene locations to 
# bridge the two datasets 
# AIM: use gene order from MGLO to set order of progressive cactus blocks

library(tidyverse)
#install.packages("comprehenr")
library(comprehenr)

#' Gene order supplied by MGLO in fasta-style format and re-formated to one
#' line prior to reading into R (3 columns: genome, chr, gene_list)
#' Here split the 3rd column (with the list of gene numbers) into a vector
#' and return it named with "genome_chr"
get_gene_order <- function(x) {
  genome = x[1]
  chr = x[2]
  listname <- paste(genome, chr, sep = "_")
  gene_list = x[3]
  gene_list = x
  gene_order_vec <- gene_list %>% 
    str_trim( side = c("both")) %>% 
    strsplit( ",") %>%
    as.vector()
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

# check number of genes in each genome
gene_list_lengths <- matrix(to_vec(for(x in tip_gene_lists) c(x[1], x[2], length(x[[3]]))), ncol = 3, byrow = TRUE)
colnames(gene_list_lengths) <- c("genome", "chr", "num_genes")
gene_list_lengths <- as_tibble(gene_list_lengths)
gene_list_lengths$chr <- as.numeric(gene_list_lengths$chr)
gene_list_lengths$num_genes <- as.numeric(gene_list_lengths$num_genes)
genes_per_genome <- gene_list_lengths %>% group_by(genome) %>% summarise(total_genes = sum(num_genes))
genes_per_genome

# read in gff files
Eam708 <- read.delim("data/Eam708_Epichloe_amarillans_NFE708_38224169_v2.gff", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
Ebr1 <- read.delim("data/Epichloe_bromicola_Nfe1.gff3", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
Ecl1605 <- read.delim("data/Epichloe_clarkii_Hl.gff3", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
Eel728 <- read.delim("data/Epichloe_elymi_728.gff3", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
Eel732 <- read.delim("data/Eel732_Epichloe_elymi_NFE732_33820330_v2.gff", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
EfeFl1 <- read.delim("data/EfFl1_v3.1.gff", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")

# subset to gene-data only
Eam708 <- Eam708[Eam708$V3 == "gene", ]
Ebr1 <- Ebr1[Ebr1$V3 == "gene",]
Ecl1605 <- Ecl1605[Ecl1605$V3 == "gene",]
Eel728 <- Eel728[Eel728$V3 == "gene",]
Eel732 <- Eel732[Eel732$V3 == "gene",]
EfeFl1 <- EfeFl1[EfeFl1$V3  == "gene",]




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
head(submaf, 20)



# check that list of genes from gff matches number of genes from maf (tip genome file)

gene_list_lengths <- function(){
  
}

# add list of genes positions to each block

# check number of blocks without genes





##### get ancestral gene lists and order blocks

# read in gene order files of ancestral genomes
anc_gen <- read.delim("data/geneorder_fmtR.out", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(anc_gen) <- c("genome", "chr", "gene_order")
# create a list of vectors of gene order named by genome and chr ("genome_chr").
anc_gene_lists <- apply(tip_gen, 1, get_gene_order)
str(gene_lists)

# order blocks for each maf 'ancestral' genome by ancestral gene list














