#!usr/bin/bash

# make CHROnicle dir structure
mkdir Epichloe Epichloe/00rawGenome Epichloe/01Genomes

# Eama708 Epichloe amarillans NFe708
perl ./scripts/def_from_gff.pl ../data/Epichloe_Transect_Genomes/E.amarillans_NFe708/annotation/Epichloe_amarillans.gff3 Eama708
perl ./scripts/make_prt.pl ../data/Epichloe_Transect_Genomes/E.amarillans_NFe708/annotation/Epichloe_amarillans.proteins.fa Eama708.def Eama708
for i in {1..7}; do echo -ne "$i\t" >> Eama708.ch; done; echo >> Eama708.ch
for i in {1..7}; do echo -ne "0\t" >> Eama708.ch; done; echo >> Eama708.ch
for i in {1..7}; do echo -ne "0\t" >> Eama708.ch; done; echo >> Eama708.ch
# sort new def file by chromosome and position
sort -k3,3 -k4,4 Eama708.rdef > Eama708.srdef
mv Eama708.srdef Eama708.def
rm Eama708.rdef
mv Eama708* Epichloe/01Genomes/

# Ebro1 
#perl ./scripts/def_from_gff.pl ../data/Epichloe_Transect_Genomes/E.bromicola_NFe1/Epichloe_bromicola_Nfe1.gff3 Ebro1
#perl ./scripts/make_prt.pl ../data/Epichloe_Transect_Genomes/E.bromicola_NFe1/Ebr1_Epichloe_bromicola_NFe1_46201037_v1.fna Ebro1.def Ebro1
#mv Ebro1* Epichloe/01Genomes/Ebro1

# Eclar1605_22
#perl ./scripts/def_from_gff.pl ../data/Epichloe_Transect_Genomes/E.clarkii_1605_22/annotation/Epichloe_clarkii_Hl.gff3 Eclar1605_22
#perl ./scripts/make_prt.pl ../data/Epichloe_Transect_Genomes/E.clarkii_1605_22/Ecl1605_22_Epichloe_clarkii_1605_22_45646793_v1.fna Eclar1605_22.def Eclar1605_22
#mv Eclar1605_22* Epichloe/01Genomes/Eclar1605_22
#data/Epichloe_Transect_Genomes/E.clarkii_1605_22/annotation/Ec_AT_rich.gff

# Eely728
perl ./scripts/def_from_gff.pl ../data/Epichloe_Transect_Genomes/E.elymi_NFe728/annotation/Epichloe_elymi_728.gff3 Eely728
perl ./scripts/make_prt.pl ../data/Epichloe_Transect_Genomes/E.elymi_NFe728/annotation/Epichloe_elymi_728.proteins.fa Eely728.def Eely728
# has mitochondria - numbered as chr 8 here
for i in {1..8}; do echo -ne "$i\t" >> Eely728.ch; done; echo >> Eely728.ch
for i in {1..8}; do echo -ne "0\t" >> Eely728.ch; done; echo >> Eely728.ch
for i in {1..8}; do echo -ne "0\t" >> Eely728.ch; done; echo >> Eely728.ch
# change chr 'm' to chr 8 in files too
sed -i 's/\tm\t/\t8\t/g' Eely728.def
sed -i 's/\tm\t/\t8\t/g' Eely728.prt
# sort new def file by chromosome and position
sort -k3,3 -k4,4 Eely728.rdef > Eely728.srdef
mv Eely728.srdef Eely728.def
rm Eely728.rdef
mv Eely728* Epichloe/01Genomes/
#data/Epichloe_Transect_Genomes/E.elymi_NFe728/annotation/Eel728_AT_rich.gff

# Eely732
#perl ./scripts/def_from_gff.pl ../data/Epichloe_Transect_Genomes/E.elymi_NFe732/Eel732_Epichloe_elymi_NFE732_33820330_v2.gff Eely732
#perl ./scripts/make_prt.pl ../data/Epichloe_Transect_Genomes/E.elymi_NFe732/annotation/Epichloe_elymi.proteins.fa  Eely732.def Eely732
#mv Eely732* Epichloe/01Genomes/
#data/Epichloe_Transect_Genomes/E.elymi_NFe732/annotation/Eel732_AT_rich.gff
#data/Epichloe_Transect_Genomes/E.elymi_NFe732/annotation/Epichloe_elymi.gff3

# Efes437
#perl ./scripts/def_from_gff.pl ../data/Epichloe_Transect_Genomes/E.festucae_E437/annotation/Epichloe_festucae_E437.gff3
#perl ./scripts/make_prt.pl ../data/Epichloe_Transect_Genomes/E.festucae_E437/
#mv Efes437* Epichloe/01Genomes/
#data/Epichloe_Transect_Genomes/E.festucae_E437/annotation/AT_rich.gff
#data/Epichloe_Transect_Genomes/E.festucae_E437/annotation/Epichloe_festucae_E437_repeats.gff

# EfesFl1 (note the only protein file available matches the v2 gff3 file)
perl ./scripts/def_from_gff.pl ../data/Epichloe_Transect_Genomes/E.festucae_Fl1/Old/EfeFl1_v2.gff3 EfesFl1
perl ./scripts/make_prt.pl ../data/Epichloe_Transect_Genomes/E.festucae_Fl1/EfFl1_proteins.faa EfesFl1.def EfesFl1

for i in {1..7}; do echo -ne "$i\t" >> EfesFl1.ch; done; echo >> EfesFl1.ch
for i in {1..7}; do echo -ne "0\t" >> EfesFl1.ch; done; echo >> EfesFl1.ch
for i in {1..7}; do echo -ne "0\t" >> EfesFl1.ch; done; echo >> EfesFl1.ch
# sort new def file by chromosome and position
sort -k3,3 -k4,4 EfesFl1.rdef > EfesFl1.srdef
mv EfesFl1.srdef EfesFl1.def
rm EfesFl1.rdef
mv EfesFl1* Epichloe/01Genomes/ 

#data/Epichloe_Transect_Genomes/E.festucae_Fl1/EfeFl1_AT_rich.gff  
#data/Epichloe_Transect_Genomes/E.festucae_Fl1/EfeFl1_Epichloe_festucae_Fl1_35023690_v2.fna
#data/Epichloe_Transect_Genomes/E.festucae_Fl1/EfFl1_transcripts.fna

# Etyp1756
#mkdir Epichloe/01Genomes/Etyp1756
#perl ./scripts/make_prt.pl ../data/Epichloe_Transect_Genomes/E.typhina_1756/Ety1756_Epichloe_typhina_1756_33870766_v3.fna
#mv Etyp1756* Epichloe/01Genomes/Etyp1756


