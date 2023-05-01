#!/bin/bash

PATH_TO_DIR=${1}

echo Configuring file system ... $(date +"%T")
mkdir ${PATH_TO_DIR}/
mkdir ${PATH_TO_DIR}/raw_fastqs
mkdir ${PATH_TO_DIR}/for_dada2
mkdir ${PATH_TO_DIR}/final_data
mkdir ${PATH_TO_DIR}/final_data/csv_output
mkdir ${PATH_TO_DIR}/final_data/rdata_output
mkdir ${PATH_TO_DIR}/final_data/logs
mkdir ${PATH_TO_DIR}/analysis_output
mkdir ${PATH_TO_DIR}/cutadapt_reports
echo Finished creating file system ... $(date +"%T")
echo Config finished! $(date +"%T")
