
## Merge MGLO and Progressive Cactus data


### DATA:
+ cactus.maf : contains tree in first lines   
+ geneorder.out  : contains MGLO output
+ gene_order.tree : MGLO tree


### AIM:
+ read in cactus.maf blocks (each paragraph starting with 'a' is a multiple sequence alignment)
+ dump blocks less than 100bp
+ identify genes in each block and pair with MGLO gene number
+ see if blocks can be ordered by MGLO gene order


