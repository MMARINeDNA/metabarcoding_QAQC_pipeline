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
quarto render Report_MURI_Module3.qmd --to html
cd ..
mv ./scripts/phyloseq_final.Rdata ./scripts/${RUN_NAME}_phyloseq_final.Rdata
mv ./scripts/${RUN_NAME}_phyloseq_final.Rdata ./final_data/rdata_output/
mv ./scripts/*html ./analysis_output/
echo finished step 3. $(date +"%T")
sleep 2
echo metabarcoding pipeline complete! $(date +"%T")

###############################################################