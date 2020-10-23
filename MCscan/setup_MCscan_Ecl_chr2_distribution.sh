#!usr/bin/bash

## rename chrs to match between gff and genome.fna files
#sed -i 's/Ecl_1605_22/Ecl1605/g' data/Epichloe_clarkii_Hl.gff3
#sed -i 's/chr/EfeFl1/g' data/EfFl1_v3.1.gff
# note formatted all fna file to be 60 chars for ease of use
# formatted all chr to match between gff and fna files

#######################################################
#            MAKE CACTUS ROOT BED AND CDS FILES       #
#######################################################
#echo "creating bed and cds files.............."

# create bedfiles from gff files
#python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Eam708_Epichloe_amarillans_NFE702_38224169_v2.gff -o data/Eam.bed
#python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_bromicola_Nfe1.gff3 -o data/Ebr.bed
# may need to put in E fes E437 
#python -m jcvi.formats.gff bed --type=mRNA --key=Name data/EfFl1_v3.1.gff -o data/Efl1.bed
#python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_elymi_728.gff3 -o data/Eel728.bed
#python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Eel732_Epichloe_elymi_NFE732_33820330_v2.gff -o data/Eel732.bed
#python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_clarkii_Hl.gff3 -o data/Ecl.bed
#python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_typhina_Dg.gff3 -o data/Ety.bed
#python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_typhina_poae_NFe76.gff -o data/Etp.bed

# extract sequences from genome fasta using bed co-ordinates
# reduce to max length transcript and then to chr:from-to
# for faidx, region file must be in the format chr:from-to
#for i in Eam Ebr Efl1 Eel728 Eel732 Ecl Ety Etp;
#do
#    Rscript --vanilla  scripts/longest_transcript.R data/$i.bed
#done

# replace old bedfile with longest transcript bedfile
#for i in  Eam Ebr Efl1 Eel728 Eel732 Ecl Ety Etp
#do
#    mv data/$i.maxlen.bed data/$i.bed
#done

# get cds data for each file using samtools
#samtools faidx data/Eam708_Epichloe_amarillans_NFE708_38224169_v2.fna  --region-file data/Eam.regions.txt  --output data/Eam.cds  --mark-strand sign
#samtools faidx data/Ebr1_Epichloe_bromicola_NFe1_46201037_v1.fna --region-file data/Ebr.regions.txt  --output data/Ebr.cds  --mark-strand sign
#samtools faidx data/EfeFl1_Epichloe_festucae_Fl1_35023690_v2.fna --region-file data/Efl1.regions.txt --output data/Efl1.cds --mark-strand sign
#samtools faidx data/Eel728_Epichloe_elymi_NFe728_34206040_v2.fna --region-file data/Eel728.regions.txt  --output data/Eel728.cds  --mark-strand sign
#samtools faidx data/Eel732_Epichloe_elymi_NFE732_33820330_v2.fna --region-file data/Eel732.regions.txt  --output data/Eel732.cds  --mark-strand sign
#samtools faidx data/Ecl1605_22_Epichloe_clarkii_1605_22_45646793_v1.fna --region-file data/Ecl.regions.txt --output data/Ecl.cds  --mark-strand sign
#samtools faidx data/Ety1756_Epichloe_typhina_1756_33870766_v3.fna --region-file data/Ety.regions.txt  --output data/Ety.cds  --mark-strand sign
#samtools faidx data/Etp76_Epichloe_typhina_var_poae_NFe76_38327242_v1.fna --region-file data/Etp.regions.txt --output data/Etp.cds --mark-strand sign

# note [faidx] Truncated sequence: Ecl1605_3:7592542-7597818 for E clarkii

# pop mRNA_# data back into cds headers from lookup file
#for i in Eam Ebr Efl1 Eel728 Eel732 Ecl Ety Etp; 
#do
#    # linearise cds file
#    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < data/$i.cds > data/$i.lin.cds
#    # replace with region lookup file info
#    sed 's/mrna/>mrna/g' data/$i.regions_lookup.txt > data/$i.regions_headers.txt
#    cut -f 2 data/$i.lin.cds > data/$i.seqs.cds
#    paste data/$i.regions_headers.txt data/$i.seqs.cds > data/$i.cds
#    # use perl to change the second instance of tab to a carrige return in each line
#    perl -pe 's{\t}{++$n % 2 ? $& : "\n"}ge' data/$i.cds > temp
#    mv temp data/$i.cds
#    # re-format to 60 chars per line
#    cp data/$i.cds .
#    R --slave --vanilla < scripts/format_60char_fasta.R --args $i.cds
#    mv fmt_$i.cds data/$i.cds
#    rm $i.cds
#    # MCscan takes the first number in a chromosome name as the chromosome number. 
#    # change chromosome names as appropriate to accommodate this
#done


########################################################
#              USE MCscan
#########################################################

#echo "create and populate results directory ................."
#mkdir results_chr2
cd results_chr2
#for i in Eam Ebr Efl1 Eel728 Eel732 Ecl Ety Etp
#do
#    cp ../data/$i.bed .
#    cp ../data/$i.cds .
#done

# format CDS files (to transcript name only)
#for i in Eam Ebr Efl1 Eel728 Eel732 Ecl Ety Etp
#do
#    python -m jcvi.formats.fasta format $i.cds > tmp.cds
#    mv tmp.cds $i.cds
#done

# synteny blocks
#echo "calculate synteny blocks ................."
#rm Eam.Ebr*
#rm Ebr.Efl1*
#rm Efl1.Eel728*
#rm Eel728.Eel732*
#rm Eel732.Ecl*
#rm Ecl.Ety*
#rm Ety.Etp*

python -m jcvi.compara.catalog ortholog Eam Ebr --no_strip_names
#python -m jcvi.compara.catalog ortholog Ebr Efl1 --no_strip_names
#python -m jcvi.compara.catalog ortholog Efl1 Eel728 --no_strip_names
#python -m jcvi.compara.catalog ortholog Eel728 Eel732 --no_strip_names
#python -m jcvi.compara.catalog ortholog Eel732 Ecl --no_strip_names
#python -m jcvi.compara.catalog ortholog Ecl Ety --no_strip_names
#python -m jcvi.compara.catalog ortholog Ety Etp --no_strip_names


# filter synteny blocks
#echo "calculate filtered synteny blocks............."
#rm .last.filtered



#python3 -m jcvi.compara.catalog ortholog Eam Ebr --cscore=.99 --no_strip_name
#python3 -m jcvi.compara.catalog ortholog Ebr Efl1 --cscore=.99 --no_strip_name
#python3 -m jcvi.compara.catalog ortholog Efl1 Eel728 --cscore=.99 --no_strip_name
#python3 -m jcvi.compara.catalog ortholog Eel728 Eel732 --cscore=.99 --no_strip_name
#python3 -m jcvi.compara.catalog ortholog Eel732 Ecl --cscore=.99 --no_strip_name
#python3 -m jcvi.compara.catalog ortholog Ecl Ety --cscore=.99 --no_strip_name
#python3 -m jcvi.compara.catalog ortholog Ety Etp --cscore=.99 --no_strip_name

# compare synteny reciprocal depth of coverage <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Doesn't work for this data set?????
#echo "compare synteny reciprocal depth of coverage.............."
#python -m jcvi.compara.synteny depth --histogram Eam.Ebr.anchors
#python -m jcvi.compara.synteny depth --histogram Ebr.Efl1.anchors
#python -m jcvi.compara.synteny depth --histogram Efl1.Eel728.anchors
#python -m jcvi.compara.synteny depth --histogram Eel728.Eel732.anchors
#python -m jcvi.compara.synteny depth --histogram Eel732.Ecl.anchors
#python -m jcvi.compara.synteny depth --histogram Ecl.Ety.anchors
#python -m jcvi.compara.synteny depth --histogram Ety.Etp.anchors

# make dotplot (already done by filter synteny block command)
#echo "make dotplot ............"
#python -m jcvi.graphics.dotplot Eam.Ebr.anchors
#python -m jcvi.graphics.dotplot Ebr.Efl1.anchors
#python -m jcvi.graphics.dotplot Efl1.Eel728.anchors
#python -m jcvi.graphics.dotplot Eel728.Eel732.anchors
#python -m jcvi.graphics.dotplot Eel732.Ecl.anchors
#python -m jcvi.graphics.dotplot Ecl.Ety.anchors
#python -m jcvi.graphics.dotplot Ety.Etp.anchors

##### manually make seqids and layout plots  ############

### seqids < list of comma separated chr names, one spp per line
#chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19
#Pp01,Pp02,Pp03,Pp04,Pp05,Pp06,Pp07,Pp08

# list chrs that map to E. clarkii chr 2 (from test_MLGO_gene_distributions)
echo '''Eam708_3,Eam708_5
Ebr1_3
EfeFl1_1,EfeFl1_2,EfeFl1_6
Eel728_2
Eel732_1
Ecl1605_2
Ety1756_2
Etp76_1,Etp76_2,Etp76_3,Etp76_4,Etp76_5,Etp76_6'''> seqids

### layout < example graph layout
## y, xstart, xend, rotation, color, label, va,  bed
# .6,     .1,    .8,       0,      , Grape, top, grape.bed
# .4,     .1,    .8,       0,      , Peach, top, peach.bed
## edges
#e, 0, 1, grape.peach.anchors.simple

echo '''# y, xstart, xend, rotation, color, label, va,  bed
.8,     .1,    .8,       0,       , Eam, top, Eam.bed
.7,     .1,    .8,       0,       , Ebr, top, Ebr.bed
.6,     .1,    .8,       0,       , Efl1, top, Efl1.bed
.5,     .1,    .8,       0,       , Eel728, top, Eel728.bed
.4,     .1,    .8,       0,       , Eel732, top, Eel732.bed
.3,     .1,    .8,       0,       , Ecl, top, Ecl.bed
.2,     .1,    .8,       0,       , Ety, top, Ety.bed
.1,     .1,    .8,       0,       , Etp, top, Etp.bed
# edges
e, 0, 1, Eam.Ebr.anchors.simple
e, 1, 2, Ebr.Efl1.anchors.simple
e, 2, 3, Efl1.Eel728.anchors.simple
e, 3, 4, Eel728.Eel732.anchors.simple
e, 4, 5, Eel732.Ecl.anchors.simple
e, 5, 6, Ecl.Ety.anchors.simple
e, 6, 7, Ety.Etp.anchors.simple''' > layout


# make karyotype plot
##echo "make karyotype plot..............."
#python3 -m jcvi.compara.synteny screen --minspan=30 --simple Ecl.Efl1.anchors Ecl.Efl1.anchors.new
#python3 -m jcvi.compara.synteny screen --minspan=30 --simple Efl1.Etp.anchors Efl1.Etp.anchors.new
#python3 -m jcvi.graphics.karyotype seqids layout

cd ..


####### SORTED figures by hand
#Note: removed Chr18,Chr87,Chr185,Chr121,Chr130,Chr197,Chr175,Chr213 from ssid list for MLGO chr 5 as they were duplicated from MLGO chr 1
 
