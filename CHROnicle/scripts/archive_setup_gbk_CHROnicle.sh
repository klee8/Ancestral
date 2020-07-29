# unzip Epichloe genbank files from David
cd GB_files_from_DW
tar -xvf epichloe_GB.tar.gz

# rename files
mv Epichloe_amarillans.gbk EamaNfE708.gbk
mv Epichloe_bromicola_Nfe1.gbk EbroNfE1.gbk
mv Epichloe_clarkii_Hl.gbk Eclar_1605_22.gbk
mv Epichloe_elymi_728.gbk EelyNfE728.gbk
mv Epichloe_elymi.gbk EelyNfE732.gbk
mv Epichloe_festucae_E437.gbk EfesE437.gbk
mv Epichloe_festucae_Fl1.gbk EfesFl1.gbk
mv Epichloe_typhina_Dg.gbk Etyp_1756.gbk
mv Epichloe_typhina_poae_NFe76.gbk EtypNfE76.gbk
cd ..


# setup Chronicle file structure for Epichloe genomes
mkdir Epichloe Epichloe/00rawGenom Epichloe/01Genomes
for i in `cat genome_list.txt`; do mkdir Epichloe/00rawGenom/$i; cp GB_files_from_DW/$i.gbk Epichloe/00rawGenom/$i/; done

# add info to ACCESSION line in gbk files so ConvertGBK.py can interact with it
for i in `cat genome_list.txt`
do 
    echo "setting up $i.dat file for ConvertGbk.py" 
    cd Epichloe/00rawGenom/$i
    # make a .dat copy of the .gbk file
    cp $i.gbk $i.dat
    # grab reference line with chr base positions and add to sed command
    grep 'REFERENCE' $i.gbk > tmp
    cut -d ' ' -f 7,9 tmp > bases
    sed -i 's/)//g' bases
    sed -i 's/ /:/' bases
    sed -i "s/$/:1\/\' $i.dat/" bases
    # create chr number list
    for j in {1..8}; do echo "/ACCESSION chromosome:KLMU1:$j:" >> chr; done
    #sed -i "s/^/\/ACCESSION chromosome:KLMU1:/" chr
    # get line numbers for accessions and add into sed command
    grep -n 'ACCESSION' $i.dat > accessions
    sed -i 's/ //g' accessions
    sed -i 's/:/s\//g' accessions
    sed -i "s/^/sed -ie \'/g" accessions
    # paste it all together
    paste --delimiters="\0" accessions chr bases > add_accession_info.sh
    # remove intermediate files
    rm accessions bases chr tmp
    bash add_accession_info.sh
    cd ../../../
done


# format input files as per instructions in CHRONInicle SynChro
cd Programs/0Convert2InputF
for i in `cat ../../genome_list.txt`; do echo "converting $i"; ./ConvertGbk.py Epichloe $i $i; done 
cd ../../
