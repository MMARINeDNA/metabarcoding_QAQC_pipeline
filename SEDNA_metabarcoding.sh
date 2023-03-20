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
# This wrapper takes in _______ as inputs and has an interactive interface that allows you to run the entire pipeline and decide what primers/databases to use for the analysis.
# USAGE: bash metabarcoding_wrapper.sh {INPUT_FASTQ}

module load R

# NOTE : NEED TO CHANGE PATHWAYS TO PATHWAYS ON SEDNA. CURRENTLY SET TO DIRECTORY STRUCTURE ON CEG
INPUT_FASTQ = ${1}
ALL_TAX_REF = ${2}
ALL_PRIMER_DATA = ~/aim/muri_project/mod3/example_dir_structure/metadata/MURI301.csv


#################### STEP 0: cd into pipeline directory and move input ####################
cd ~/aim/muri_project/mod3/example_dir_structure
mv ./fastq_holding_pen/${INPUT_FASTQ} ./raw_fastqs/
###############################################################

#################### STEP 1: DADA2 ####################

read -p 'Taxonomy Database to Input?:  ' first_tax
RScript ./scripts/dada2QAQC.R ./raw_fastqs/${INPUT_FASTQ} ${first_tax} ${ALL_PRIMER_DATA} 
mv ./raw_fastqs/*.Rdata ../for_more_tax

###############################################################

#################### STEP 2: Assign Taxonomy (again) ####################
while true; do
read -p 'Another round of Taxonomy? (y/n) ' tax_reply
if [[ ${tax_reply} == "y" ]]; then
cat assignTaxonomy_codes
read -p "Enter Primer Name: " primer

# mifish
if [[ $primer = "MFU" ]]
then
RScript ./scripts/assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}


# marver 1
elif [[ $primer = "MV1" ]]
then
RScript ./scripts/assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}


# d-loop
elif [[ $primer = "DL" ]]
then
RScript ./scripts/assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}


# C18 
elif [[ $primer = "C18" ]]
then
RScript ./scripts/assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}

# custom database input
else 
echo "Seems like that isn't a primer we recognize..."
read -p 'Input a custom database: ' custom_db
if [ -z ${custom_db} ]
break
else 
RScript ./scripts/assignTaxonomy.R ./for_more_tax/*.Rdata ${custom_db}
fi
fi

else
break
fi
done

mv 

###############################################################

#################### STEP 3: Generate Report ####################

Rscript -e "rmarkdown::render('Primer_test_prelim_report.qmd',params=list(args = myarg))"

###############################################################
