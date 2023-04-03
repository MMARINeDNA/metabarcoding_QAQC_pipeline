#!/bin/bash

#SBATCH --partition=node # Queue selection
#SBATCH --job-name=build_reference_database_sedna# Job name
#SBATCH --mail-type=ALL # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user= [EMAIL HERE] # Where to send mail\
#SBATCH --ntasks=1 # Run a single task\
#SBATCH --cpus-per-task=20  # Number of CPU cores per task\
#SBATCH --mem=10000 # Job memory request\
#SBATCH --time=24:00:00 # Time limit hrs:min:sec\
#SBATCH --output=build_reference_database.log # Standard output/error\
#export OMP_NUM_THREADS=8


# Metabarcoding Wrapper
# written by: Alexandria Im
# This wrapper takes in run you want to run as an input and has an interactive interface that allows you to run the entire pipeline and decide what primers/databases to use for the analysis.
# USAGE: bash metabarcoding_wrapper.sh {FIRST_TAX} {SECOND_TAX}

# NOTE : NEED TO CHANGE PATHWAYS TO PATHWAYS ON SEDNA. CURRENTLY SET TO DIRECTORY STRUCTURE ON CEG

#################### VARIABLE ASSIGNMENT ####################
FIRST_TAX=${1}
SECOND_TAX=${2}
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
rm ./raw_fastqs/*
echo finished primer and quality trimming. $(date +"%T")
sleep 3

###############################################################

#################### STEP 2: DADA2 ####################

echo starting step 1: dada2 ... $(date +"%T")
RScript ./scripts/dada2QAQC.R ./for_dada2 ${FIRST_TAX} ${ALL_PRIMER_DATA} 
mv ./*.Rdata ./for_more_tax
rm ./for_dada2/*
echo finished step 1. $(date +"%T")
sleep 3

###############################################################

#################### STEP 3: Assign Taxonomy (again) ####################

echo starting step 2: assigning taxonomy again... $(date +"%T")
sleep 2
if [ -z ${SECOND_TAX} ]; then
break

# mifish
elif [[ ${SECOND_TAX} = "MFU" ]]
then
RScript ./scripts/assignTaxonomy.R ./metadata/MiFish_12S_0223_dada2.fasta ${primer}


# marver 1
elif [[ ${SECOND_TAX} = "MV1" ]]
then
RScript ./scripts/assignTaxonomy.R ./metadata/MiFish_12S_0223_dada2.fasta ${primer}


# d-loop
elif [[ ${SECOND_TAX} = "DL" ]]
then
RScript ./scripts/assignTaxonomy.R ./metadata/cetacean_DL_taxonomy.fasta ${primer}


# C16 
elif [[ ${SECOND_TAX} = "C16" ]]
then
RScript ./scripts/assignTaxonomy.R ./metadata/ceph_C16_sanger.fasta ${primer}
fi


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







