#!/bin/bash

# Metabarcoding Wrapper
# written by: Alexandria Im
# This wrapper takes in run you want to run as an input and has an interactive interface that allows you to run the entire pipeline and decide what primers/databases to use for the analysis.
# USAGE: bash metabarcoding_wrapper.sh {INPUT_FASTQ}

# NOTE : NEED TO CHANGE PATHWAYS TO PATHWAYS ON SEDNA. CURRENTLY SET TO DIRECTORY STRUCTURE ON CEG

#################### VARIABLE ASSIGNMENT ####################
ALL_PRIMER_DATA=~/Desktop/muri_sandbox/example_data_structure/metadata/MURI301.csv
MFU_F="GCCGGTAAAACTCGTGCCAGC"
MFU_R="CATAGTGGGGTATCTAATCCCAGTTTG"
DL_F="TCACCCAAAGCTGRARTTCTA"
DL_R="GCGGGTTGCTGGTTTCACG"
MV1_F="CGTGCCAGCCACCGCG"
MV1_R="GGGTATCTAATCCYAGTTTG"
C16_F="GACGAGAAGACCCTAWTGAGCT"
C16_R="AAATTACGCTGTTATCCCT"

#################### STEP 0: cd into pipeline directory and move input ####################
cd ~/Desktop/muri_sandbox/example_data_structure/
# mv ./fastq_holding_pen/${INPUT_FASTQ} ./raw_fastqs/
###############################################################

#################### STEP 1: Use Cutadapt ####################
echo starting primer and quality trimming... $(date +"%T")
sleep 3
cd raw_fastqs
CUTADAPT=$(which cutadapt)
for i in *R1*
do
FILE_PRIM=$(echo ${i} | cut -d - -f 1) #grab primer name at beginning of file name
FILE_NAME=$(echo ${i} | cut -d _ -f 1,2,3) #grab everything before R1 (aka sample name)
R1=${i}
R2=$(echo ${i} | sed 's/R1/R2/g')
if [[ ${FILE_PRIM} == "MFU" ]]; then
echo MFU Detected
${CUTADAPT} -g ${MFU_F} \
     -G "${MFU_R}" \
     -o ../for_dada2/${R1} \
     -p ../for_dada2/${R2} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
"${R1}" "${R2}" 1> "../cutadapt_reports/${FILE_NAME}_trim_report.txt"
elif [[ ${FILE_PRIM} == "DL" ]]; then
echo DL Detected
${CUTADAPT} -g ${DL_F} \
     -G "${DL_R}" \
     -o ../for_dada2/${R1} \
     -p ../for_dada2/${R2} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
"${R1}" "${R2}" 1> "../cutadapt_reports/${FILE_NAME}_trim_report.txt"
elif [[ ${FILE_PRIM} == "MV1" ]]; then
echo MV1 Detected
${CUTADAPT} -g ${MV1_F} \
     -G "${MV1_R}" \
     -o ../for_dada2/${R1} \
     -p ../for_dada2/${R2} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
"${R1}" "${R2}" 1> "../cutadapt_reports/${FILE_NAME}_trim_report.txt"
elif [[ ${FILE_PRIM} == "C16" ]]; then
echo MV1 Detected
${CUTADAPT} -g ${C16_F} \
     -G "${C16_R}" \
     -o ../for_dada2/${R1} \
     -p ../for_dada2/${R2} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
"${R1}" "${R2}" 1> "../cutadapt_reports/${FILE_NAME}_trim_report.txt"

fi
done

#grabbing important info from cutadapt reports and synthesize in overall_report.txt
echo starting reports...
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
echo finished primer and quality trimming. $(date +"%T")
sleep 3

###############################################################

#################### STEP 2: DADA2 ####################

echo starting step 1: dada2 ... $(date +"%T")
sleep 3
read -p 'Taxonomy Database to Input?:  ' first_tax
RScript ./scripts/dada2QAQC.R ./for_dada2 ${first_tax} ${ALL_PRIMER_DATA} 
mv ./*.Rdata ./for_more_tax
echo finished step 1. $(date +"%T")
sleep 3

###############################################################

#################### STEP 3: Assign Taxonomy (again) ####################

echo starting step 2: assigning taxonomy again... $(date +"%T")
sleep 2
while true; do
read -p 'Another round of Taxonomy? (y/n) ' tax_reply
if [[ ${tax_reply} == "y" ]]; then
cat ./metadata/assignTaxonomy_codes
read -p "Enter Ref DB Name: " primer

# mifish
if [[ $primer = "MFU" ]]
then
RScript ./scripts/assignTaxonomy.R ./metadata/MiFish_12S_0223_dada2.fasta ${primer}


# marver 1
elif [[ $primer = "MV1" ]]
then
RScript ./scripts/assignTaxonomy.R ./metadata/MiFish_12S_0223_dada2.fasta ${primer}


# d-loop
elif [[ $primer = "DL" ]]
then
RScript ./scripts/assignTaxonomy.R ./metadata/cetacean_DL_taxonomy.fasta ${primer}


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

echo finished step 2. $(date +"%T")
sleep 3

###############################################################

#################### STEP 4: Generate Report ####################
echo starting step 3: making the stats file... $(date +"%T")
sleep 3
cd ./scripts
Rscript -e "rmarkdown::render('Primer_test_prelim_report.qmd')"
cd ..
mv ./scripts/primer_test_phyloseq.Rdata ./final_data/
mv ./scripts/*html ./analysis_output/
echo finished step 3. $(date +"%T")
sleep 2
echo metabarcoding pipeline complete! $(date +"%T")

###############################################################

# TRIALS AND TRIBULATIONS
# what does truncQ do? how different from cutadapt? what is the scale for scoring (2 is so low???)
# is it cutadapt or dada2 that can't handle I's?
# 

# TO-DO:
# change output for assignTax directory IN R CODE
# do we want it to ask all questions at the beginning of the run? or is it ok to ask when step 2 starts for example?
# asv/hashing for dada2 step
# holding pen will be on QNAP and rest will be on SEDNA... need to figure out how to automatically move
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
# change for loop in dada2 script to not break when it hits a primer that isn't there
# enabling inputs for the R scripts

# THINGS TO CHANGE BEFORE HAND-OFF:
# change script hard-coded paths







