#!/bin/bash

# Metabarcoding Wrapper
# written by: Alexandria Im
# This wrapper takes in _______ as inputs and has an interactive interface that allows you to run the entire pipeline and decide what primers/databases to use for the analysis.
# USAGE: bash metabarcoding_wrapper.sh {INPUT_FASTQ}

# NOTE : NEED TO CHANGE PATHWAYS TO PATHWAYS ON SEDNA. CURRENTLY SET TO DIRECTORY STRUCTURE ON CEG
# ALL_TAX_REF = ${1}
# ALL_PRIMER_DATA=~/Desktop/muri_sandbox/example_data_structure/metadata/MURI301.csv


#################### STEP 0: cd into pipeline directory and move input ####################
cd ~/Desktop/muri_sandbox/example_data_structure/
# mv ./fastq_holding_pen/${INPUT_FASTQ} ./raw_fastqs/
###############################################################

#################### STEP 1: DADA2 ####################

read -p 'Taxonomy Database to Input?:  ' first_tax
RScript ./scripts/dada2QAQC.R ./raw_fastqs ${first_tax} ${ALL_PRIMER_DATA} 
mv ./*.Rdata ./for_more_tax

###############################################################

#################### STEP 2: Assign Taxonomy (again) ####################
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

###############################################################

#################### STEP 3: Generate Report ####################
cd ./scripts
Rscript -e "rmarkdown::render('Primer_test_prelim_report.qmd',params=list(args = myarg))"
cd ..
mv ./scripts/primer_test_phyloseq.Rdata ./final_data/
mv ./scripts/*html ./analysis_output/

###############################################################


# TO-DO:
# BIG ISSUE: DOES PROMPTS WORK WHEN YOU SUBMIT SCRIPT AS A SLURM JOB OH NO
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
# for cutadapt, make first step and make sure I's are turned to N's in csv file
# really gotta change the trimming on DADA2- it's either trimming way too little (for the longer ones), and there's A LOT of varaition run-to-run

# DONE:
# gotta work on WHERE the outputs go and the general directory structure
# ask again and again to do step 2 until you're done
# ask if we want to do the assign taxonomy step again
# ask for multiple inputs

# THINGS TO CHANGE BEFORE HAND-OFF:
# change script hard-coded paths







