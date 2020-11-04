#!usr/bin/bash

## rename chrs to match between gff and genome.fna files
sed -i 's/Ecl_1605_22/Ecl1605/g' data/Epichloe_clarkii_Hl.gff3
sed -i 's/chr/EfeFl1/g' data/EfFl1_v3.1.gff
# NOTES:
# < formatted all fna file to be 60 chars for ease of use
# < formatted all chr to match between gff and fna files
# < MCscan takes the first number in a chromosome name as the chromosome number.
#   (change chromosome names as appropriate to accommodate this)


#######################################################
#            MAKE CACTUS ROOT BED AND CDS FILES       #
#######################################################
#echo "creating bed and cds files.............."

# create bedfiles from gff files
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Eam708_Epichloe_amarillans_NFE702_38224169_v2.gff -o data/Eam.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_bromicola_Nfe1.gff3 -o data/Ebr.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_festucae_E437.gff3 -o data/Efe437.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/EfFl1_v3.1.gff -o data/Efl1.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_elymi_728.gff3 -o data/Eel728.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Eel732_Epichloe_elymi_NFE732_33820330_v2.gff -o data/Eel732.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_clarkii_Hl.gff3 -o data/Ecl.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_typhina_Dg.gff3 -o data/Ety.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name data/Epichloe_typhina_poae_NFe76.gff -o data/Etp.bed

# extract sequences from genome fasta using bed co-ordinates
# reduce to max length transcript and then to chr:from-to
# longest_transcript.R outputs *.regions.txt files for samtools faidx
# it also outputs *.regions.lookup.txt files to format cds headers later
# for faidx, region file must be in the format chr:from-to
for i in Eam Ebr Efe437 Efl1 Eel728 Eel732 Ecl Ety Etp
do
    Rscript --vanilla  scripts/longest_transcript.R data/$i.bed
done

# replace old bedfile with longest transcript bedfile
for i in  Eam Ebr Efe437 Efl1 Eel728 Eel732 Ecl Ety Etp
do
    mv data/$i.maxlen.bed data/$i.bed
done

# get cds data for each file using samtools
samtools faidx data/Eam708_Epichloe_amarillans_NFE708_38224169_v2.fna  --region-file data/Eam.regions.txt  --output data/Eam.cds  --mark-strand sign
samtools faidx data/Ebr1_Epichloe_bromicola_NFe1_46201037_v1.fna --region-file data/Ebr.regions.txt  --output data/Ebr.cds  --mark-strand sign
samtools faidx data/EfeFl1_Epichloe_festucae_Fl1_35023690_v2.fna --region-file data/Efl1.regions.txt --output data/Efl1.cds --mark-strand sign
samtools faidx data/EfeE437_Epichloe_festucae_E437_33219473_v1.fna --region-file data/Efe437.regions.txt --output data/Efe437.cds --mark-strand sign
samtools faidx data/Eel728_Epichloe_elymi_NFe728_34206040_v2.fna --region-file data/Eel728.regions.txt  --output data/Eel728.cds  --mark-strand sign
samtools faidx data/Eel732_Epichloe_elymi_NFE732_33820330_v2.fna --region-file data/Eel732.regions.txt  --output data/Eel732.cds  --mark-strand sign
samtools faidx data/Ecl1605_22_Epichloe_clarkii_1605_22_45646793_v1.fna --region-file data/Ecl.regions.txt --output data/Ecl.cds  --mark-strand sign
samtools faidx data/Ety1756_Epichloe_typhina_1756_33870766_v3.fna --region-file data/Ety.regions.txt  --output data/Ety.cds  --mark-strand sign
samtools faidx data/Etp76_Epichloe_typhina_var_poae_NFe76_38327242_v1.fna --region-file data/Etp.regions.txt --output data/Etp.cds --mark-strand sign

# note [faidx] Truncated sequence: Ecl1605_3:7592542-7597818 for E clarkii

# pop mRNA_# data back into cds headers from lookup file
for i in Eam Ebr Efe437 Efl1 Eel728 Eel732 Ecl Ety Etp; 
do
    # linearise cds file
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < data/$i.cds > data/$i.lin.cds
    # replace with region lookup file info
    sed 's/mrna/>mrna/g' data/$i.regions_lookup.txt > data/$i.regions_headers.txt
    cut -f 2 data/$i.lin.cds > data/$i.seqs.cds
    paste data/$i.regions_headers.txt data/$i.seqs.cds > data/$i.cds
    # use perl to change the second instance of tab to a carrige return in each line
    perl -pe 's{\t}{++$n % 2 ? $& : "\n"}ge' data/$i.cds > temp
    mv temp data/$i.cds
    # re-format to 60 chars per line
    cp data/$i.cds .
    R --slave --vanilla < scripts/format_60char_fasta.R --args $i.cds
    mv fmt_$i.cds data/$i.cds
    rm $i.cds
    # MCscan takes the first number in a chromosome name as the chromosome number. 
    # change chromosome names as appropriate to accommodate this
done

