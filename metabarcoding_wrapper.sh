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
sh ../scripts/primer_trimming.sh

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
RScript ./scripts/dada2QAQC.R ./for_dada2/ ./final_data/ ./metadata/ ${RUN_NAME}
echo finished step 1. $(date +"%T")
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
# is it cutadapt lsor dada2 that can't handle I's?
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







