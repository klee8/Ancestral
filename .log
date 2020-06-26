  

## AIM: Reconstruct ancestral genomes

Testing output of different programs to see if we can reconstruct fungal genome ancestor sensibly

### DATA

MLGO output of gene order
Progressive Cactus from alignment data

#### Cactus data
DATA from Amelie
Export the HAL alignment as MAF referenced at the root
```
hal2maf cactus.hal cactus.maf
```
#### MGLO data
DATA from Murray Cox:
+ ortholog files
  1 The gene order in modern species is in ‘gene_order.txt’ 
  2	the phylogeny was fixed from standard phylogeny methods ‘gene_order.phy’
+ MGLO files
- geneorder.out
- gene_order.dat
geneorder.out file has lists of gene orders for ancestral genomes. 
ancestral genome A6 is the genome of interest

### Checked the distribution of E.clarkii chromosome 2 genes on Ancestral genome 6 (root of tree)
Expect that it should be contiguous on the ancestor based on distribution at the tips
First check how many of the genes on each E.clarkii chromosomes end up on each ancestral chromosome
```
perl check_mglo_output_overlaps.pl Ecl_gene_order.dat A6_geneorder.out > Ecl_chr_distribution_on_A6_genome.txt
```
Graph data to see how they compare across tip genomes
```
R slave --vanilla < 

```

### Checked the contiguous blocks of E.clarkii chromosome 2 genes on all tip and ancestral genomes
Looked at the presence/absence (1/0) of E.clarkii chr 2 genes on each tip and ancestral genome.
(Note that this did not include checking the gene order, just pres/abs).

```
# tips
perl flag_gene_blocks.pl Ecl_chr2.txt Ecl_gene_order.dat MGLO_tagged.out
# ancestors
perl flag_gene_blocks.pl Ecl_chr2.txt geneorder.out MGLO_tagged.out
```








