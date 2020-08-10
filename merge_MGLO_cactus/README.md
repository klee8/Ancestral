
## Merge MGLO and Progressive Cactus data

### DATA:
+ cactus.maf : contains tree in first lines   
+ gene_order.dat : MGLO numbered genes in each tip genome
+ geneorder.out  : MGLO gene number order for each ancestral chr
+ gene_order.tree : MGLO tree
+ gff : gff files for each tip genome

### AIM:
+ read in cactus.maf blocks (each paragraph starting with 'a' is a multiple sequence alignment)
+ dump blocks less than 100bp
+ read in gene_order.dat
+ read in gff files for each tip genome
+ pair up gene order numbers with positional info and gene names for each species
+ identify genes in each cactus block and pair with MGLO gene number
+ see if blocks can be ordered by MGLO gene order in ancestral geneome (A6 in cactus, A0 in MGLO).




