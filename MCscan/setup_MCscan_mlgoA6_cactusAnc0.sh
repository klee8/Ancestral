###!usr/bin/bash
#!/bin/sh

#######################################################
#            MAKE CACTUS ROOT BED AND CDS FILES       #
#######################################################
echo "creating cactus root ancestral genome bed and cds files.............."

# create bedfile for Eam from gff file 
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Eam708_Epichloe_amarillans_NFE702_38224169_v2.gff -o data/Eam.bed

# create bedfile for Anc0 (cactus root ancestor) from Eam bed file
halLiftover data/cactus.hal Eam data/Eam.bed Anc0 data/cactus_Anc0.bed

# extract Anc0 fasta file from hal file
hal2fasta data/cactus.hal Anc0 > data/cactus_Anc0.fasta
hal2fasta data/cactus.hal Eam > data/Eama702.fasta

# extract sequences from Anc0 fasta using bed co-ordinates
# reduce to max length transcript and then to chr:from-to
# for faidx, region file must be in the format chr:from-to
Rscript --vanilla scripts/longest_transcript.R data/Eam.bed
samtools faidx data/cactus_Anc0.fasta --region-file data/cactus_Anc0_regions.txt --output data/cactus_Anc0.cds --mark-strand sign
# save the original bed file, clean up new bed file and replace old one
mv data/cactus_Anc0.bed data/cactus_Anc0_all.bed
mv data/cactus_Anc0.maxlen.bed data/cactus_Anc0.bed
sed -i '/V1 V2/d' data/cactus_Anc0.bed

# pop mRNA_# data back into cds headers from lookup file
# linearise cds file
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < data/cactus_Anc0.cds > data/cactus_Anc0.lin.cds
# replace with region lookup file info
sed 's/mrna/>mrna/g' data/cactus_Anc0_regions_lookup.txt > data/cactus_Anc0_regions_headers.txt
cut -f 2 data/cactus_Anc0.lin.cds > data/cactus_Anc0.seqs.cds
paste data/cactus_Anc0_regions_headers.txt data/cactus_Anc0.seqs.cds > data/cactus_Anc0.cds
# use perl to change the second instance of tab to a carrige return in each line
perl -pe 's{\t}{++$n % 2 ? $& : "\n"}ge' data/cactus_Anc0.cds > temp
mv temp data/cactus_Anc0.cds
# re-format to 60 chars per line
cp data/cactus_Anc0.cds .
R --slave --vanilla < scripts/format_60char_fasta.R --args cactus_Anc0.cds
mv fmt_cactus_Anc0.cds data/cactus_Anc0.cds
rm cactus_Anc0.cds
# MCscan takes the first number in a chromosome name as the chromosome number. 
# remove "Anc0ref" from bed and cds files to prevent all chromosomes being labelled as '0'
sed -i 's/Anc0ref//g' data/cactus_Anc0.cds
sed -i 's/Anc0ref//g' data/cactus_Anc0.bed



#######################################################
#            MAKE MLGO ROOT BED AND CDS FILES       #
#######################################################
echo "creating MLGO root ancestral genome bed and cds files................."

# extract Eam fasta file from hal file
hal2fasta data/cactus.hal Eam > data/Eama702.fasta

# create MLGO ancestral genome (A6) bedfile and list regions and headers for fasta
R --slave --vanilla < scripts/setup_MLGO_bed_fasta_files.R 
sed -i '/gen_chr/d' data/MLGO_A6.bed

# setup MLGO fasta
sed -i '/^x/d' data/MLGO_A6_Eama_regions.txt
samtools faidx data/Eama702.fasta
samtools faidx data/Eama702.fasta --region-file data/MLGO_A6_Eama_regions.txt --output data/MLGO_A6.cds --mark-strand sign
# linearise cds file
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < data/MLGO_A6.cds > data/MLGO_A6.lin.cds
# replace with region lookup file info
sed -i '/^x/d' data/MLGO_A6_fasta_headers.txt
cut -f 2 data/MLGO_A6.lin.cds > data/MLGO_A6.seqs.cds
paste data/MLGO_A6_fasta_headers.txt data/MLGO_A6.seqs.cds > data/MLGO_A6.cds
sed -i 's/\t/\n/g' data/MLGO_A6.cds
# re-format to 60 chars per line
cp data/MLGO_A6.cds .
R --slave --vanilla < scripts/format_60char_fasta.R --args MLGO_A6.cds
mv fmt_MLGO_A6.cds data/MLGO_A6.cds
rm MLGO_A6.cds



########################################################
#              USE MCscan
#########################################################

echo "create and populate results directory ................."
mkdir results
cd results
cp ../data/cactus_Anc0.bed .
cp ../data/cactus_Anc0.cds .
cp ../data/MLGO_A6.bed .
cp ../data/MLGO_A6.cds .

# make bed files
#python -m jcvi.formats.gff bed --type=mRNA --key=Name Vvinifera_145_Genoscope.12X.gene.gff3.gz -o grape.bed
#python -m jcvi.formats.gff bed --type=mRNA --key=Name Ppersica_298_v2.1.gene.gff3.gz -o peach.bed

# format CDS files (to transcript name only)
#python -m jcvi.formats.fasta format Vvinifera_145_Genoscope.12X.cds.fa.gz grape.cds
#python -m jcvi.formats.fasta format Ppersica_298_v2.1.cds.fa.gz peach.cds

# synteny blocks
#echo "calculate synteny blocks ................."
#python -m jcvi.compara.catalog ortholog MLGO_A6 cactus_Anc0 --no_strip_names

# filter synteny blocks
echo "calculate filtered synteny blocks............."
#rm MLGO_A6.cactus_Anc0.last.filtered
python3 -m jcvi.compara.catalog ortholog MLGO_A6 cactus_Anc0 --cscore=.99 --no_strip_name

# compare synteny reciprocal depth of coverage <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Doesn't work for this data set?????
#echo "compare synteny reciprocal depth of coverage.............."
#python -m jcvi.compara.synteny depth --histogram MLGO_A6.cactus_Anc0.anchors

# make dotplot (already done by filter synteny block command)
#echo "make dotplot ............"
#python -m jcvi.graphics.dotplot MLGO_A6.cactus_Anc0.anchors


##### manually make seqids and layout plots  ############

### seqids < list of comma separated chr names, one spp per line
#chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19
#Pp01,Pp02,Pp03,Pp04,Pp05,Pp06,Pp07,Pp08

cut -f 1 MLGO_A6.bed | sort | uniq > tmp
for i in `cat tmp`; do echo -ne "$i," >> seqids; done
echo >> seqids
cut -f 1 cactus_Anc0.bed | sort | uniq > tmp
# to sort numerically temporarily remove 'Chr' from chr name
sed 's/Chr//g' tmp | sort -n | sed 's/^/Chr/g' > temp
for i in `cat temp`; do echo -ne "$i," >> seqids; done
rm tmp temp
sed -i 's/,$//g' seqids 

### layout < example graph layout
## y, xstart, xend, rotation, color, label, va,  bed
# .6,     .1,    .8,       0,      , Grape, top, grape.bed
# .4,     .1,    .8,       0,      , Peach, top, peach.bed
## edges
#e, 0, 1, grape.peach.anchors.simple

echo '''# y, xstart, xend, rotation, color, label, va,  bed
.6,     .1,    .8,       0,       , MLGO, top, MLGO_A6.bed
.4,     .1,    .8,       0,       , cactus, top, cactus_Anc0.bed
# edges
e, 0, 1, MLGO_A6.cactus_Anc0.anchors.simple''' > layout


# make karyotype plot
echo "make karyotype plot..............."
python3 -m jcvi.compara.synteny screen --minspan=30 --simple MLGO_A6.cactus_Anc0.anchors MLGO_A6.cactus_Anc0.anchors.new
python3 -m jcvi.graphics.karyotype seqids layout

cd ..


####### SORTED figures by hand
#Note: removed Chr18,Chr87,Chr185,Chr121,Chr130,Chr197,Chr175,Chr213 from ssid list for MLGO chr 5 as they were duplicated from MLGO chr 1
 
