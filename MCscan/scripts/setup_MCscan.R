# setup MCscan

library(tidyverse)
#library(stringr)

########################################################
#   Create MLGO 'bed' and 'fasta' files for MCscan     #
########################################################

# note that MLGO does not have positional information for genes
# take ortholog numbers from MLGO, map them to Eama with flat ortho file
# use Eama mRNA sequence position from gff use these to get mRNA 
# nucleotide sequences and put a fake 1KB ATATAT sequence between
# each CDS for the MLGO 'fasta' file for each chromosome

#' Gene order supplied by MGLO in fasta-style format and re-formated to one
#' line prior to reading into R (3 columns: genome, chr, gene_list)
#' Here split the 3rd column (with the list of gene numbers) into a vector
#' and return it as a dataframe with AncGen_chr, orthogroup and strand columns
split_mlgo_to_lists <- function(x) {
  genome = x[1]
  chr = x[2]
  listname <- paste(genome, chr, sep = "_")
  gene_list = x
  gene_order_vec <- gene_list %>% 
    str_trim( side = c("both")) %>% 
    strsplit( ",")
  return(gene_order_vec)
}

#' convert nested list of orthogroup gene names to dataframe
#' add in column with ancestral gene and chr name pasted together
#' add in column of strand info
#' strip strand info from orthogroupname
nested_list_to_dataframe <- function(nestedlist) {
    df <- data.frame(matrix(unlist(nestedlist[[1]][[3]]), nrow=length(nestedlist[[1]][[3]]), byrow=T))
    colnames(df) <- c("orthogroup")
    df$gen_chr <- rep(paste(nestedlist[[1]][[1]], nestedlist[[1]][[2]], sep = "_"), length(nestedlist[[1]][[3]]))
    df$strand <- ifelse(grepl("-",df$orthogroup), "-", "+")
    df$orthogroup <- gsub("-", "", df$orthogroup)
  return(df)
}

# map Eama genes to ortho numbers in mlgo ancestral genomes
mlgo_AncGen <- read.delim("data/geneorder.fmtR.out", header = FALSE, stringsAsFactors = FALSE)
colnames(mlgo_AncGen) <- c("AncGen", "chr", "orthlist")
orthos <- read.delim("data/ortho_long.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(orthos) <- c("orthogroup", "transcript", "species")

#' take in MGLO data for one chromosome
#' split into a list of "genome", "chr", "list of geneorders"
#' re-arrange into a dataframe ("gen_chr", "strand", "orthogroup")
mlgo_to_df <- function(mlgo_AncGen, genome, chromosome) {
  temp <- mlgo_AncGen[mlgo_AncGen$AncGen == genome & mlgo_AncGen$chr == chromosome, ]   
  temp2 <- apply(temp, 1, split_mlgo_to_lists)
  df <- nested_list_to_dataframe(temp2)
  return(df)
}

# for each mlgo chromosome gene list for A6 ancient genome
# create a dataframe with "orthogroup", "genome_chr" and "strand"
# bind dataframes for each chromosome together for whole genome 
df <- NULL    # setup an empty object to append to 
for(chromosome in 1:6){
  tempdf <- mlgo_to_df(mlgo_AncGen, "A6", chromosome)
  df <- rbind(df, tempdf)  
}
#view(df)

# preserve mglo ordering of orthogroups 
ordering <- structure(.Names = df$orthogroup, 1:length(df$orthogroup))

# merge mglo dataframe with Eam702 info from flat ortho file
orthos <- orthos[orthos$species == "Eam702",]
orthos$orthogroup <- gsub("og_0*", "", orthos$orthogroup, perl = TRUE)
df <- merge(df, orthos, by = "orthogroup")

# merge Eam702 gene positions from gff to df <<<<<<<<<<<<<< (note Eam strands do not match the ancestral A6 ones)
Eam702gff <- read.delim("data/Eam708_Epichloe_amarillans_NFE702_38224169_v2.gff", header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
Eam702gff <- Eam702gff[Eam702gff$V3 == "mRNA",]
Eam702gff$V9 <- sub(";.*", "", sub("ID=", "", Eam702gff$V9), perl = TRUE)
df <- merge(df, Eam702gff, by.x = "transcript", by.y = "V9")

# re-order dataframe by mlgo ancestral genome order
df <- df[order(ordering[df$orthogroup]), ]

#' add in fake positions for MLGO orthologs 
#' (mapped as Eam genes) and set 1kb apart
add_fake_mlgo_gene_positions_end <- function(df) {
  df <- df %>% group_by(gen_chr) %>% mutate(acc_length = cumsum(length))                
  df$inc <- rep(1000, length(df$length))
  df <- df %>% group_by(gen_chr) %>% mutate(acc_inc = cumsum(inc))
  df$end <- c(1, rep(0, length(df$length)-1))
  df$end <- ifelse(!(df$end == 1), (lag(df$acc_inc) + df$acc_length), df$length)
  df <- subset(df, select = -c(inc, acc_inc, acc_length) )
  return(df)
}

#' Add in fake start positions
#' assumes end positions are already calculated
add_fake_mlgo_gene_positions_start<- function(df) {
  df$start <- c(1, rep(0, length(df$end)-1))
  df$start <- ifelse(!(df$start == 1), lag(df$end)+1000, 1)   
  return(df)
}

# add bedfile elements                  
#add_fake_bed_type <- function(df){                                        #<<<<<<<<<<<<<< are you supposed to reset mRNA numbers for new chrs >>> NOPE!!!
#  df$type <- paste("mRNA", 1:length(df$length), sep = "_")
#  return(df)
#}

# calculate the length of each mRNA
df$length <- abs(df$V5-df$V4)

# for each chromosome, add in fake MLGO ortholog positions and bedfile type
# making sure to re-set start/end/mRNA number for each line
tempdf2 <- NULL   # setup an empty object to append to 
for(chromosome in 1:6){
  tempdf <- df[df$gen_chr == paste("A6", chromosome, sep = "_"), ]
  tempdf <- add_fake_mlgo_gene_positions_end(tempdf)
  tempdf <- add_fake_mlgo_gene_positions_start(tempdf)
  tempdf2 <- bind_rows(tempdf2, tempdf)  
}
df <- tempdf2

# add in bed file mRNA name (just numbered in sequence) and score (between 0-1000) 
df$name <- paste("mRNA", 1:length(df$length), sep = "_")
df$score <- rep(0, length(df$name))

# subset to bedfile
bedfile <- df[, c("gen_chr", "start", "end", "name", "score", "strand")]

# write to file
write.table(bedfile, "MGLO_A6.bed", quote = FALSE, row.names = FALSE)
 
