#!/bin/bash

# Metabarcoding Wrapper
# written by: Alexandria Im
# This wrapper takes in _______ as inputs and has an interactive interface that allows you to run the entire pipeline and decide what primers/databases to use for the analysis.
# USAGE: bash metabarcoding_wrapper.sh {INPUT_FASTQ}

# NOTE : NEED TO CHANGE PATHWAYS TO PATHWAYS ON SEDNA. CURRENTLY SET TO DIRECTORY STRUCTURE ON CEG
INPUT_FASTQ = ${1}
ALL_TAX_REF = ~/aim/muri_project/mod3/example_dir_structure/metadata/NCBI_mtDNA_alltaxa_reference.fasta
ALL_PRIMER_DATA = ~/aim/muri_project/mod3/example_dir_structure/metadata/MURI301.csv
TAX_REF_2 = ${4}

#################### STEP 0: cd into pipeline directory and move input ####################
cd ~/aim/muri_project/mod3/example_dir_structure
mv ./fastq_holding_pen/${INPUT_FASTQ} ./raw_fastqs/
###############################################################

#################### STEP 1: DADA2 ####################

RScript dada2QAQC.R ./raw_fastqs/${INPUT_FASTQ} ${ALL_TAX_REF} ${ALL_PRIMER_DATA} 
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
RScript assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}


# marver 1
elif [[ $primer = "MV1" ]]
then
RScript assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}


# marver 3
elif [[ $primer = "MV3" ]]
then
RScript assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}


# d-loop
elif [[ $primer = "DL" ]]
then
RScript assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}


# leray
elif [[ $primer = "LRY" ]]
then
RScript assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}


# C16
elif [[ $primer = "C16" ]]
then
RScript assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}


# C18 
elif [[ $primer = "C18" ]]
then
RScript assignTaxonomy.R ./for_more_tax/*.Rdata ${TAX_REF_2}

else 
echo idk that one... did you make a typo?
fi

else
break
fi
done

mv 

###############################################################

#################### STEP 3: Generate Report ####################

Rscript -e "rmarkdown::render('example.Rmd',params=list(args = myarg))"

###############################################################


# TO-DO:
# look into if you can do rmd in this context?
# do we want it to ask all questions at the beginning of the run? or is it ok to ask when step 2 starts for example?
# asv/hashing for dada2 step
# enabling inputs for the R scripts

# look into why we wanna 
# for trimming, might be trimming way too much for the longer ones like Dloop (429 trimmed to 230)
# look into calculating where to trim 

# DONE:
# gotta work on WHERE the outputs go and the general directory structure
# ask again and again to do step 2 until you're done
# ask if we want to do the assign taxonomy step again
# ask for multiple inputs








