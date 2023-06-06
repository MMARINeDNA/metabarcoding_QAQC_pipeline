#!/bin/bash

PATH_TO_DIR=${1}

echo Configuring file system ... $(date +"%T")
mkdir ${PATH_TO_DIR}/raw_fastqs
mkdir ${PATH_TO_DIR}/for_dada2
mkdir ${PATH_TO_DIR}/final_data
mkdir ${PATH_TO_DIR}/final_data/csv_output
mkdir ${PATH_TO_DIR}/final_data/rdata_output
mkdir ${PATH_TO_DIR}/final_data/logs
mkdir ${PATH_TO_DIR}/analysis_output
mkdir ${PATH_TO_DIR}/cutadapt_reports


echo Copying over necessary files ... $(date +"%T")
cp  -r ./metadata/* ${PATH_TO_DIR}/metadata
cp   ./bin/* ${PATH_TO_DIR}/scripts
echo Config finished! $(date +"%T")
