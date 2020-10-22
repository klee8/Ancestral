# merge Progressive Cactus and MGLO output using gff gene locations to 
# bridge the two datasets 
# AIM: use gene order from MGLO to set order of progressive cactus blocks

library(tidyverse)
#install.packages("comprehenr")
library(comprehenr)



############################################################
###    Get MLGO ordered gene lists for extant species   ####
############################################################

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

# read in gene order files of tip genomes
tip_gen <- read.delim("data/gene_order_fmtR.dat", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(tip_gen) <- c("genome", "chr", "gene_order")
tipnames <- paste(tip_gen$genome, tip_gen$chr, sep="_")

# create a list of vectors of gene order named by genome and chr ("genome_chr").
tip_gene_lists <- apply(tip_gen, 1, get_gene_order)
names(tip_gene_lists) <- tipnames

# check number of genes in each genome
gene_list_lengths <- matrix(to_vec(for(x in tip_gene_lists) c(x[1], x[2], length(x[[3]]))), ncol = 3, byrow = TRUE)
colnames(gene_list_lengths) <- c("genome", "chr", "num_genes")
gene_list_lengths <- as_tibble(gene_list_lengths)
gene_list_lengths$chr <- as.numeric(gene_list_lengths$chr)
gene_list_lengths$num_genes <- as.numeric(gene_list_lengths$num_genes)
genes_per_genome <- gene_list_lengths %>% group_by(genome) %>% summarise(total_genes = sum(num_genes))
genes_per_genome

#genome  total_genes
#1 Eam708         7451
#2 Ebr1           7615
#3 Ecl1605        7772
#4 Eel728         7423
#5 Eel732         8323
#6 EfeFl1         7822
#7 Etp76          7504
#8 Ety1756        7539


###########################################################
###  GFF FILES: find gene positions in extant species   ###
###########################################################

#' subset and merge gff information
#' 
#' takes in gff file
#' subsets to gene data only
#' reduces to columns of gene name, chromosome, start, end and strand
#' appends new gff to master gff file
subset_and_merge_gff <- function(new_gff, master_gff) {
  new_gff <- new_gff[new_gff$V3 == "gene",]
  new_gff <- new_gff[,c("V9", "V1", "V4", "V5", "V7")]
  new_gff$V9 <- sub("ID=", "", new_gff$V9)
  new_gff$V9 <- sub(";", "", new_gff$V9)
  colnames(new_gff) <- c("gene_names", "chr", "start_pos", "end_pos", "strand")
  master_gff <- rbind(master_gff, new_gff)
  return(master_gff)
}

# read in gff files
Eam708 <- read.delim("data/Eam708_Epichloe_amarillans_NFE702_38224169_v2.gff", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
Ebr1 <- read.delim("data/Epichloe_bromicola_Nfe1.gff3", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
Ecl1605 <- read.delim("data/Epichloe_clarkii_Hl.gff3", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
Eel728 <- read.delim("data/Epichloe_elymi_728.gff3", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
Eel732 <- read.delim("data/Eel732_Epichloe_elymi_NFE732_33820330_v2.gff", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
EfeFl1 <- read.delim("data/EfeFl1_v3.gff", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
Etp76 <- read.delim("data/Epichloe_typhina_poae_NFe76.gff", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
Ety1756 <- read.delim("data/Epichloe_typhina_Dg.gff3", stringsAsFactors = FALSE, header = FALSE, comment.char = "#")

master_gff <- data.frame(gene_names=character(), start_pos=integer(), end_pos=integer(), strand=factor() )
head(Eam708)
master_gff <- subset_and_merge_gff(Eam708, master_gff)
master_gff <- subset_and_merge_gff(Ebr1, master_gff)
master_gff <- subset_and_merge_gff(Ecl1605, master_gff)
master_gff <- subset_and_merge_gff(Eel728, master_gff)
master_gff <- subset_and_merge_gff(Eel732, master_gff)
master_gff <- subset_and_merge_gff(EfeFl1, master_gff)
master_gff <- subset_and_merge_gff(Etp76, master_gff)
master_gff <- subset_and_merge_gff(Ety1756, master_gff)
str(master_gff)

# check that list of genes from gff matches number of genes from maf (tip genome file)  <<< they don't
#nrow(Eam708)    # 7451
#nrow(Ebr1)      # 7738
#nrow(Ecl1605)   # 7881
#nrow(Eel728)    # 7541
#nrow(Eel732)    # 8324
#nrow(EfeFl1)    # 7926
#nrow(Etp76)     # 7620
#nrow(Ety1756)   # 7660


###########################################################################################
###  Use ortholog flat file to map gene positions from gff to ortholog numbers in MLGO  ###
###########################################################################################

# read in ortho_long.tsv file (where the numbers came from)
orthos <- read.delim("data/ortho_long.tsv", stringsAsFactors = FALSE, header = FALSE)
colnames(orthos) <- c("ortholog", "transcript", "species")
head(orthos)

# merge with 'numeric' gene names from MGLO
orthos$transcript <- sub("-T1", "", orthos$transcript)
orth_position <- merge(orthos, master_gff, by.x="transcript", by.y="gene_names")
orth_position$ortholog <- sub("og_", "", orth_position$ortholog)
str(orth_position)


#######################################################################
###                      Get cactus maf file                        ###
#######################################################################

# read in re-formatted cactus maf file  (from `perl fmt_cactus_for_R.pl cactus.maf`)
maf <- read.delim("data/cactus.fmtR.maf", stringsAsFactors = FALSE, header = FALSE)
colnames(maf) <- c("block", "s", "genome.chr", "pos", "length", "strand", "chrlength", "sequence")
# remove s column
maf$s <- NULL
submaf <- maf %>% separate(genome.chr, into = c("genome", "chr"), sep = "\\.")

#### Look to see properties of maf file
## total blocks <<< 424422
#dim(Eam_maf)
## total length of blocks <<< 21306977
#sum(Eam_maf$length)
## number of blocks with genes  <<<< 6110
#dim(Eam_maf[Eam_maf$orthlist != "none",])
## number of blocks without genes <<<< 418312
#dim(Eam_maf[Eam_maf$orthlist == "none",])
## summed length of blocks with genes   <<<< 2427254
#sum(Eam_maf[Eam_maf$orthlist != "none", c("length")])
## summed length of blocks without genes   <<<< 18879723
#sum(Eam_maf[Eam_maf$orthlist == "none", c("length")])
#Eam_maf

# check the numbers for blocks of length >= 500
#sub_Eam_maf <- Eam_maf[Eam_maf$length >= 500, ]    
#dim(sub_Eam_maf)            # 9,514
#sum(sub_Eam_maf$length)     # 7,817,467
#dim(sub_Eam_maf[sub_Eam_maf$orthlist != "none",])      # 2,051
#sum(sub_Eam_maf[sub_Eam_maf$orthlist != "none", c("length")]) #  1,741,557
#unique(Eam_maf$chr)


#####   GRAPHICAL summary of Eam blocks

sum_Eam_blocks <- as_tibble(Eam_maf) %>% 
  mutate( 
    genehere = as.factor(ifelse(orthlist != "none", 1, 0)), 
    genenothere = ifelse(orthlist == "none", 1, 0),
    lengenehere = ifelse(orthlist != "none", length, 0),
    lengenenothere = ifelse(orthlist == "none", length, 0)) 

head(sum_Eam_blocks)
ggplot(data = sum_Eam_blocks) +
  geom_col(mapping = aes(x = chr, y = length, fill = genehere))
ggsave("Total_length_of_blocks_with_genes_per_chr.pdf")


# check that list of genes in each block does not clash with order in genome order vectors from MGLO
str(tip_gene_lists)     #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# reduce to Eam lines only
Eam_maf <- submaf %>% filter(genome == regex("Eam")) %>%
  select(block, genome, chr, pos, length, strand, chrlength)
# take first instance of Eam on each block (can have two or more, especially for shorter blocks)
Eam_maf <- Eam_maf[!duplicated(Eam_maf$block),]


####################################################################################
###    Identify genes in cactus maf file blocks based on Eam block positions    ####
####################################################################################

#' return true if vector is of length zero
vector_is_empty <- function(x) return(length(x) ==0 )

#' Take in vector of orthologs and return concatenated string or "none"
#'  
#' genelist <- as.list(c("780", "900"))
#' list_of_genes_to_string(genelist)
#' > "708,900"
#' genelist <- as.list(c())
#' list_of_genes_to_string(genelist)
#' > "none"
list_of_genes_to_string <- function(genelist){
  genestring <- { if(vector_is_empty(genelist)) "none" else paste(genelist, collapse = ",") }
  genestring <- { if(genestring == "character(0)") "none" else paste(genelist, collapse = ",") }
  return(genestring)
}

#' Take in block attributes (chromosome, position and length)
#' Take in positions of orthologs 'orth_position'
#' return list of orthologs that have a start position in the block
select_genes_in_block <- function(chromosome, position, length) {
  genelist <- orth_position %>% 
    filter(chr == chromosome & start_pos >= position & start_pos <= (position + length)) %>% 
    select(ortholog)
  genestring <- list_of_genes_to_string(genelist, sep = ",")
  return(genestring)
}

# test mapply output
#temp <- mapply(select_genes_in_block, head(Eam_maf$chr, n=300), head(Eam_maf$pos, n=300), head(Eam_maf$length, n=300))
#temp2 <- as.tibble(temp) %>% filter(value != "none")
#view(temp2)

# map list of genes to each block based on block positions in Eam sequences
Eam_maf$orthlist <- mapply(select_genes_in_block, Eam_maf$chr, Eam_maf$pos, Eam_maf$length)


#####################################################################################
###    Subset cactus maf file to root ancestor and merge with block gene info    ####
#####################################################################################

# subset maf file to ancestral 0 genome
Anc0_maf <- submaf %>% filter(genome == regex("Anc0")) %>%
  select(block, genome, chr, pos, length, strand, chrlength)

## check lengths of Anc0 and Eam maf files (note, don't have the samenumber of lines. Anc0 can have 0-many Eam alignments)
#length(unique(Anc0_maf$block))  # 477275
#length(unique(Eam_maf$block))   # 415977

# add Eam genelists to Anc0 maf file 
AncGenOrthlist <- merge(Anc0_maf, Eam_maf[, c("block", "orthlist")], by="block", all = TRUE)

# get a list of genes on each Anc0 chromosome (where they are available)
Anc0Genes <- AncGenOrthlist %>% group_by(chr) %>%
  filter(!is.na(orthlist) & orthlist!="none") %>% 
  summarise(list_of_genes_to_string(orthlist))
colnames(Anc0Genes) <- c("AncChr", "geneList")

# split string of genes back into a list
Anc0Genes <- Anc0Genes %>% mutate(list = strsplit(geneList, split = ","))
head(Anc0Genes)


#################################################################################
###   get ancestral gene lists from MLGO and check against gene list in maf   ###
#################################################################################

# read in gene order files of ancestral genomes
anc_gen <- read.delim("data/geneorder.fmtR.out", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(anc_gen) <- c("genome", "chr", "gene_order")

# select root ancestor
MLGOroot <- anc_gen[anc_gen$genome == "A6", ]
# get rid of all strand info on gene list
MLGOroot$gene_order <- gsub("-", "", MLGOroot$gene_order)
# split gene_order string into a list
MLGOroot <- MLGOroot %>% mutate(MLGOgenelist = strsplit(gene_order, split = ","))


################################################################################
###                    WRITE OUT FILES TO PROCESS WITH PERL                  ###
################################################################################


write.table(MLGOroot[,c("genome", "chr","gene_order")],"MLGOroot.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(Anc0Genes[,c("AncChr", "geneList")], "Anc0Genes", quote = FALSE, row.names = FALSE, sep = "\t")


#################################################################################
###             Summarising perl-processed comparison file                    ###
#################################################################################

# Perl processed file took in gene lists from MGLO (root ancestor A6) and from cactus
# root ancestor (Anc0) and indexes them by chromosome and order
# it subset the MGLO list to genes that were on a particular cactus chromosome
# and then goes throught that chromosome in cactus order and gives the MLOG chromosme
# and index position
# headers: Gene   CactusChr CactusIndex MLGOchr MLGOindex

compared <- read.delim("MLGO_cactus_compared.txt", header = TRUE, stringsAsFactors = FALSE)

# for each cactus chromosome, get a list of MLGO chromosomes that map to it and the number of them
conflict <- compared %>% group_by(CactusChr) %>% summarise(numMLGOchrs =  length(unique(MLGOchr)),
                                                           MLGOchrs = paste(unique(MLGOchr), collapse = ","))

# find the number of cactus chromosomes with more than one MLGO chromsome genes on it
dim(conflict[conflict$numMLGOchrs > 1, ])

# count the number of each MLGO chromosome combination against cactus chromsomes and 
# count their occurances 
conflict %>% group_by(MLGOchrs) %>% summarise(num = n())

# look at the combinations of MLGO chromosomes that turn up on cactus chromosomes
unique(conflict[conflict$numMLGOchrs > 1, c("MLGOchrs") ])


###### GRAPH 
conflict
head(compared)

listMLGO <- compared %>% group_by(CactusChr) %>% summarise( list = gsub(",", "", list_of_genes_to_string(MLGOchr)))
listMLGO
write.table(listMLGO, "cactus_chrs_MLGO_gene_overlaps.txt", quote = FALSE, row.names = FALSE)
# suggested functions: ?rle geom_segments geom_rect







  
















########################################################################
###              Compare MLGO and cactus gene lists                  ###<<<<<<<<<<<<<  NOT WORKING YET
########################################################################

compare_lists <- function(MLGOdf, chromosome, CactusGeneList) {
  MLGOcheck <- MLGOdf[MLGOdf$chr == chromosome,] %>% filter(MLGOgenelist %in% CactusGeneList)
  ans <- ifelse(length(MLGOcheck$MLGOlist) < length(list), "extra", 
                ifelse(MLGOcheck$MLGOlist == goodlist, "same", "don't know"))
  #return(MLGOcheck)
  return(ans)
}

# TEST
#REFdf <- data.frame(MLGOlist = c(3,5,8,10,12))
#goodlist <- c(5,8,12)  # subset in the right order
#badorderlist <- c(3,12,8)  # wrong order
#extrageneslist <- c(3,4,8,10) # extra genes
#list = goodlist
#list = badorderlist
#list = extrageneslist
#MLGOcheck <- REFdf %>% filter(MLGOlist %in% list)
#ifelse(length(MLGOcheck$MLGOlist) < length(list), "extra genes", 
#       ifelse(MLGOcheck$MLGOlist == goodlist, "same order", "don't know"))

Anc0Genes$list
Anc0Genes <- mutate(Anc0Genes, MLGOchr1 = compare_lists(MLGOroot, 1, Anc0Genes$list),
                               MLGOchr2 = compare_lists(MLGOroot, 2, Anc0Genes$list),
                               MLGOchr3 = compare_lists(MLGOroot, 3, Anc0Genes$list),
                               MLGOchr4 = compare_lists(MLGOroot, 4, Anc0Genes$list),
                               MLGOchr5 = compare_lists(MLGOroot, 5, Anc0Genes$list),
                               MLGOchr6 = compare_lists(MLGOroot, 6, Anc0Genes$list))


################################################################################################
###### check that MLGO root geneome (A6) gene order matches those in maf root (Anc0)

#' Takes in reference dataframe (colnames = ("genome", "chr", "gene_order") and a list of genes
#' reduces ref dataframe to items on list only and drops extra factors
#' returns boolean list whether items are the same and in the same order
compare_lists <- function(MLGOdf, chromosome, CactusGeneList) {
  MLGOcheck <- MLGOdf[MLGOdf$chr == chromosome,] %>% filter(MLGOgenelist %in% CactusGeneList)
  #MLGOcheck <- {if(vector_is_empty(MLGOcheck)) "none" else MLGOcheck}
  #MLGOcheck <- {if(MLGOcheck == "none") "none" else droplevels(MLGOcheck) }
  str(data.frame(CactusGeneList))
  str(MLGOcheck)
  #return ({if(MLGOcheck == "none") "none" else data.frame(CactusGeneList) == MLGOcheck })
}

compare_lists(MLGOroot, 1, Anc0Genes$list[[3]])
colnames(MLGOroot)

Anc0Genes <- Anc0Genes %>% mutate(list = strsplit(geneList, split = ","))

Anc0Genes$list[[3]]
str(MLGOroot)
MLGOcheck <- MLGOroot[MLGOroot$chr == 1, ] %>% filter(MLGOgenelist %in% Anc0Genes$list[[3]] )
MLGOcheck <- droplevels(MLGOcheck)
str(MLGOcheck)




#########################################################################################
###   order blocks for each maf 'ancestral' genome by MLGO ancestral gene list       ### 
#########################################################################################




