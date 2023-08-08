#!/bin/bash

# Metabarcoding Wrapper
# written by: Alexandria Im
# This wrapper takes in run you want to run as an input and has an interactive interface that allows you to run the entire pipeline and decide what primers/databases to use for the analysis.
# USAGE: bash metabarcoding_wrapper.sh {path to muri_metabarcoding directory} {run name}

# NOTE : NEED TO CHANGE PATHWAYS TO PATHWAYS ON SEDNA. CURRENTLY SET TO DIRECTORY STRUCTURE ON CEG

#################### VARIABLE ASSIGNMENT ####################
INPUT_DIR=${1}
RUN_NAME=${2}

#################### STEP 0: cd into pipeline directory and move input ####################
cd ${INPUT_DIR}
###############################################################

#################### STEP 1: Use Cutadapt ####################
echo starting primer trimming... $(date +"%T")
sleep 3
cd raw_fastqs
bash ../scripts/primer_trimming.sh

#grabbing important info from cutadapt reports and synthesize in cutadapt_overall_report.txt
echo starting reports... $(date +"%T")
cd ../cutadapt_reports
for i in *
do 
FILE_NAME=$(echo ${i} | cut -d . -f 1)
echo ${FILE_NAME} >> cutadapt_overall_report.txt 
grep Quality-trimmed $i >> cutadapt_overall_report.txt 
grep "Reads with adapters" $i >> cutadapt_overall_report.txt 
grep "Reads written (passing filters):" $i >> cutadapt_overall_report.txt 
done

mv cutadapt_overall_report.txt ../final_data/logs/
cd ..
echo finished primer trimming. $(date +"%T")
sleep 3

###############################################################

#################### STEP 2: DADA2 ####################

echo starting step 1: dada2 ... $(date +"%T")
sleep 3
RScript ./scripts/dada2QAQC.R ./for_dada2 ./final_data ./metadata ${RUN_NAME}
echo finished step 1. $(date +"%T")
sleep 3

###############################################################

#################### STEP 4: Generate Report ####################
echo starting step 3: making the stats file... $(date +"%T")
sleep 3
list_markers=$(ls ./final_data/rdata_output | cut -d _ -f 5 | cut -d . -f 1 | uniq)
cd ./scripts
for marker in ${list_markers}
do
quarto render Report_MURI_Module3.qmd -P marker:${marker} --to html 
mv ./phyloseq_final.Rdata ../final_data/rdata_output/${RUN_NAME}_${marker}_phyloseq_final.Rdata
mv ./Report_MURI_Module3.html ../analysis_output/${RUN_NAME}_${marker}_Report_MURI_Module3.html
done
cd ..
echo finished step 3. $(date +"%T")
sleep 2
echo metabarcoding pipeline complete! $(date +"%T")

###############################################################







