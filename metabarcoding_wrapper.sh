#!/bin/bash

# Metabarcoding Wrapper
# written by: Alexandria Im
# This wrapper takes in run you want to run as an input and has an interactive interface that allows you to run the entire pipeline and decide what primers/databases to use for the analysis.
# USAGE: bash metabarcoding_wrapper.sh {INPUT_FASTQ}

# NOTE : NEED TO CHANGE PATHWAYS TO PATHWAYS ON SEDNA. CURRENTLY SET TO DIRECTORY STRUCTURE ON CEG

#################### VARIABLE ASSIGNMENT ####################
ALL_PRIMER_DATA=~/Desktop/muri_sandbox/example_data_structure/metadata/MURI301.csv
MFU_F="GTCGGTAAAACTCGTGCCAGC"
MFU_R="CATAGTGGGGTATCTAATCCCAGTTTG"
DL_F="TCACCCAAAGCTGRARTTCTA"
DL_R="GCGGGTTGCTGGTTTCACG"


#################### STEP 0: cd into pipeline directory and move input ####################
cd ~/Desktop/muri_sandbox/example_data_structure/
# mv ./fastq_holding_pen/${INPUT_FASTQ} ./raw_fastqs/
###############################################################

#################### STEP 1: Use Cutadapt ####################
cd raw_fastqs
CUTADAPT=$(which cutadapt)
for i in *
do
FILE_PRIM=$(echo ${i} | cut -d - -f 1)
FILE_NAME=$(echo ${i} | cut -d . -f 1)
if [[ ${FILE_PRIM} == "MFU" ]]; then
if echo ${i} | grep -q "R1" ; then
echo MFU R1 Detected
${CUTADAPT} -a ${MFU_F} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
${i}> ../for_dada2/${i} 2> ../cutadapt_reports/${FILE_NAME}_trim_report.txt
else
echo MFU R2 Detected
${CUTADAPT} -a ${MFU_R} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
${i}> ../for_dada2/${i} 2> ../cutadapt_reports/${FILE_NAME}_trim_report.txt
fi
elif [[ ${FILE_PRIM} == "DL" ]]; then
if echo ${i} | grep -q "R1" ; then
echo DL R1 Detected
${CUTADAPT} -a ${DL_F} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
${i}> ../for_dada2/${i} 2> ../cutadapt_reports/${FILE_NAME}_trim_report.txt
else
echo DL R2 Detected
${CUTADAPT} -a ${DL_R} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
${i}> ../for_dada2/${i} 2> ../cutadapt_reports/${FILE_NAME}_trim_report.txt
fi
fi
done

#grabbing important info from cutadapt reports and synthesize in overall_report.txt
cd ../cutadapt_reports
for i in *
do 
FILE_NAME=$(echo ${i} | cut -d . -f 1)
echo ${FILE_NAME} >> overall_report.txt 
grep Quality-trimmed $i >> overall_report.txt 
grep "Reads with adapters" $i >> overall_report.txt 
grep "Reads written (passing filters):" $i >> overall_report.txt 
done
cd ..

###############################################################

#################### STEP 2: DADA2 ####################

echo starting step 1 
read -p 'Taxonomy Database to Input?:  ' first_tax
RScript ./scripts/dada2QAQC.R ./for_dada2 ${first_tax} ${ALL_PRIMER_DATA} 
mv ./*.Rdata ./for_more_tax
echo finished step 1

###############################################################

#################### STEP 3: Assign Taxonomy (again) ####################

echo starting step 2
while true; do
read -p 'Another round of Taxonomy? (y/n) ' tax_reply
if [[ ${tax_reply} == "y" ]]; then
cat ./metadata/assignTaxonomy_codes
read -p "Enter Primer Name: " primer

# mifish
if [[ $primer = "MFU" ]]
then
RScript ./scripts/assignTaxonomy.R ./metadata/MiFish_12S_0223_dada2.fasta ${primer}


# marver 1
elif [[ $primer = "MV1" ]]
then
RScript ./scripts/assignTaxonomy.R ${TAX_REF_2} ${primer}


# d-loop
elif [[ $primer = "DL" ]]
then
RScript ./scripts/assignTaxonomy.R ./metadata/cetacean_dloop_taxonomy.fasta ${primer}


# C16 
elif [[ $primer = "C16" ]]
then
RScript ./scripts/assignTaxonomy.R ./metadata/ceph_C16_sanger.fasta ${primer}

# custom database input
else 
echo "Seems like that isn't a primer we recognize..."
read -p 'Input a custom database: ' custom_db
read -p 'Primer Name: ' ${p_name}
if [ -z ${custom_db} ]
then
break
else 
RScript ./scripts/assignTaxonomy.R ${custom_db} ${p_name}
fi
fi

else
break
fi
done

echo finished step 2

###############################################################

#################### STEP 4: Generate Report ####################
echo starting step 3: making the stats file
cd ./scripts
Rscript -e "rmarkdown::render('Primer_test_prelim_report.qmd')"
cd ..
mv ./scripts/primer_test_phyloseq.Rdata ./final_data/
mv ./scripts/*html ./analysis_output/
echo finished making the stats file

###############################################################

# TRIALS AND TRIBULATIONS
# what does truncQ do? how different from cutadapt? what is the scale for scoring (2 is so low???)
# is it cutadapt or dada2 that can't handle I's?
# 

# TO-DO:
# change output for assignTax directory IN R CODE
# do we want it to ask all questions at the beginning of the run? or is it ok to ask when step 2 starts for example?
# asv/hashing for dada2 step
# change for loop in dada2 script to not break when it hits a primer that isn't there
# enabling inputs for the R scripts
# holding pen will be on QNAP and rest will be on SEDNA... need to figure out how to automatically move
# note that DADA2
# DBS FOR MIFISH AND MARVER1 WILL BE SAME
# FINAL RUN USING CORRECT DBS AND FILE STRUCTURE
# think about hashing for DADA2 because you can then use hash key and annoted hash table to classify already-found hashes (saves time and computation)


# DONE:
# gotta work on WHERE the outputs go and the general directory structure
# ask again and again to do step 2 until you're done
# ask if we want to do the assign taxonomy step again
# ask for multiple inputs
# for cutadapt, make first step and make sure I's are turned to N's in csv file
# really gotta change the trimming on DADA2- it's either trimming way too little (for the longer ones), and there's A LOT of varaition run-to-run

# THINGS TO CHANGE BEFORE HAND-OFF:
# change script hard-coded paths







