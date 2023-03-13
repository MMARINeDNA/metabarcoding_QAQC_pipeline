#!/bin/bash
INPUT_FASTQ = ${1}
TAX_REF = ${2}
PRIMER_DATA = ${3}
OUTPUT_DIR = ${4}

TAX_REF_2 = ${5}

# STEP 1: DADA2
RScript dada2QAQC.R ${INPUT_FASTQ} ${TAX_REF} ${PRIMER_DATA} ${OUTPUT_DIR}
mv MURI_primer_test_mastertax_dada2_QAQC_outputDL.Rdata ${INPUT_FASTQ}_dada2_QAQC_outputDL.Rdata

# STEP 2: Assign Taxonomy
#ask user if want to perform this step
RScript assignTaxonomy.R ${INPUT_FASTQ}_dada2_QAQC_outputDL.Rdata ${TAX_REF_2}

# STEP 3: Generate Report


# gotta work on WHERE the outputs go and the general directory structure










