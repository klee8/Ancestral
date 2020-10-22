#!usr/bin/bash

# make CHROnicle dir structure
mkdir Epichloe Epichloe/00rawGenome Epichloe/01Genomes

function format_for_CHRONicle() {
    # sort new def file by chromosome and position
    sort -k3,3 -k4,4n $1.rdef > $1.srdef
    mv $1.srdef $1.def
    rm $1.rdef
    # add header line
    echo -e "type\tname\tchr\tstart\tend\tstrand\tsens\tIDg/chr\tIDg/all\tIDf/all" > tmp1
    cat tmp1 $1.def > tmp2
    mv tmp2 $1.def
    rm tmp1
    # sort prt file by chromosome and position
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $1.prt > temp1
    cat temp1 | sort -k 9n > temp2
    cat temp2 | sed 's/\t/\n/9' > $1.prt
    rm temp1 temp2 
    # reformat prt file with 60 char per line and a star at the end of proteins
    R --slave --vanilla < scripts/format_60char_fasta.R --args $1.prt
    mv fmt_$1.prt $1.prt
    # make ch file
    perl ./scripts/make_ch.pl $1.def $1
    mv $1* Epichloe/01Genomes/
}


# Eama708 Epichloe amarillans NFe708
#perl ./scripts/def_from_gff.pl ../data/data_from_Murray/Transect_Epichloe_Genomes/E.amarillans_NFe708/Epichloe_amarillans.gff3 Eama708
#perl ./scripts/make_prt.pl ../data/data_from_Murray/Epichloe_Transect_Genomes/E.amarillans_NFe708/annotation/Epichloe_amarillans.proteins.fa Eama708.def Eama708
perl ./scripts/def_from_gff.pl re-formatted_Epichloe_Genomes/E.amarillans_NFe708/Eam708_Epichloe_amarillans_NFE708_38224169_v2.gff EAMA
perl ./scripts/make_prt.pl re-formatted_Epichloe_Genomes/E.amarillans_NFe708/Epichloe_amarillans.proteins.fa EAMA.def EAMA
format_for_CHRONicle "EAMA"

# Ebro1 
#../data/data_from_Murray/Epichloe_Transect_Genomes/
#perl ./scripts/def_from_gff.pl ../data/data_from_Murray/Epichloe_Transect_Genomes/E.bromicola_NFe1/Epichloe_bromicola_Nfe1.gff3 EBRO
#perl ./scripts/make_prt.pl ../data/data_from_Murray/Epichloe_Transect_Genomes/E.bromicola_NFe1/Ebr1_Epichloe_bromicola_NFe1_46201037_v1.fna EBRO.def EBRO
#perl ./scripts/def_from_gff.pl re-formatted_Epichloe_Genomes/E.bromicola_NFe1/Epichloe_bromicola_Nfe1.gff3 EBRO
#perl ./scripts/make_prt.pl re-formatted_Epichloe_Genomes/E.bromicola_NFe1/Ebr1_Epichloe_bromicola_NFe1_46201037_v1.fna EBRO.def EBRO
#format_for_CHRONicle "EBRO"

# Eclar1605_22
perl ./scripts/def_from_gff.pl re-formatted_Epichloe_Genomes/E.clarkii_1605_22/Epichloe_clarkii_Hl.gff3 ECLR
perl ./scripts/make_prt.pl re-formatted_Epichloe_Genomes/E.clarkii_1605_22/Epichloe_clarkii_Hl.proteins.fa ECLR.def ECLR
format_for_CHRONicle "ECLR"

# Eely728
#perl ./scripts/def_from_gff.pl ../data/data_from_Murray/Epichloe_Transect_Genomes/E.elymi_NFe728/Epichloe_elymi_728.gff3 EEAT
#perl ./scripts/make_prt.pl ../data/data_from_Murray/Epichloe_Transect_Genomes/annotation/E.elymi_NFe728/Epichloe_elymi_728.proteins.fa EEAT.def EEAT
perl ./scripts/def_from_gff.pl re-formatted_Epichloe_Genomes/E.elymi_NFe728/Epichloe_elymi_728.gff3 EEAT
perl ./scripts/make_prt.pl re-formatted_Epichloe_Genomes/E.elymi_NFe728/Epichloe_elymi_728.proteins.fa EEAT.def EEAT
# has mitochondria - numbered as chr 8 here
## change chr 'm' to chr 8 in files too
sed -i 's/\tm\t/\t8\t/g' EEAT.def
sed -i 's/\tm\t/\t8\t/g' EEAT.prt
# one transcript is present in the fasta file but not in the def file... not sure why - remove this
sed -i '/T1/d' EEAT.prt
format_for_CHRONicle "EEAT"

# Eely732
perl ./scripts/def_from_gff.pl re-formatted_Epichloe_Genomes/E.elymi_NFe732/Eel732_Epichloe_elymi_NFE732_33820330_v2.gff EETW
perl ./scripts/make_prt.pl re-formatted_Epichloe_Genomes/E.elymi_NFe732/Epichloe_elymi.proteins.fa  EETW.def EETW
# one transcript is present in the fasta file but not in the def file... not sure why - remove this
sed -i '/T1/d' EETW.prt
format_for_CHRONicle "EETW" 

# Efes437
perl ./scripts/def_from_gff.pl re-formatted_Epichloe_Genomes/E.festucae_E437/Epichloe_festucae_E437.gff3 EFES
perl ./scripts/make_prt.pl re-formatted_Epichloe_Genomes/E.festucae_E437/Epichloe_festucae_E437.proteins.fa EFES.def EFES
format_for_CHRONicle "EFES"

# EfesFl1 (note the only protein file available matches the v2 gff3 file)
perl ./scripts/def_from_gff.pl  re-formatted_Epichloe_Genomes/E.festucae_Fl1/EfFl1_v3.1.gff EFLO
#perl ./scripts/make_prt.pl ../data/data_from_Murray/Epichloe_Transect_Genomes/E.festucae_Fl1/EfFl1_proteins.faa FFL1.def EFLO
perl ./scripts/make_prt.pl re-formatted_Epichloe_Genomes/E.festucae_Fl1/EfeFl1_proteins_v3.fa EFLO.def EFLO
format_for_CHRONicle "EFLO"

# Etyp1756
perl ./scripts/def_from_gff.pl re-formatted_Epichloe_Genomes/E.typhina_1756/Epichloe_typhina_Dg.gff3 ETYP
perl ./scripts/make_prt.pl re-formatted_Epichloe_Genomes/E.typhina_1756/Epichloe_typhina_Dg.proteins.fa ETYP.def ETYP
format_for_CHRONicle "ETYP" 

# E.typhina_var_poae_NFe76
#perl ./scripts/def_from_gff.pl re-formatted_Epichloe_Genomes/E.typhina_var_poae_NFe76/Epichloe_typhina_poae_NFe76.gff TYVP
#perl ./scripts/make_prt.pl re-formatted_Epichloe_Genomes/E.typhina_var_poae_NFe76/Epichloe_typhina_poae_NFe76.proteins.fa TYVP.def TYVP
#format_for_CHRONicle "TYVP" 

