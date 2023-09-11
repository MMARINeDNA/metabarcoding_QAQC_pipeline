#!/bin/bash

PATH_TO_DIR=${1}

echo Configuring file system ... $(date +"%T")
mkdir ${PATH_TO_DIR}/muri_metabarcoding/
mkdir ${PATH_TO_DIR}/muri_metabarcoding/raw_fastqs
mkdir ${PATH_TO_DIR}/muri_metabarcoding/for_dada2
mkdir ${PATH_TO_DIR}/muri_metabarcoding/final_data
mkdir ${PATH_TO_DIR}/muri_metabarcoding/final_data/csv_output
mkdir ${PATH_TO_DIR}/muri_metabarcoding/final_data/rdata_output
mkdir ${PATH_TO_DIR}/muri_metabarcoding/final_data/logs
mkdir ${PATH_TO_DIR}/muri_metabarcoding/final_data/logs/cutadapt_reports
mkdir ${PATH_TO_DIR}/muri_metabarcoding/analysis_output
mkdir ${PATH_TO_DIR}/muri_metabarcoding/metadata
mkdir ${PATH_TO_DIR}/muri_metabarcoding/scripts

echo Copying over necessary files ... $(date +"%T")
cp  -r ./metadata/* ${PATH_TO_DIR}/muri_metabarcoding/metadata
cp   ./bin/* ${PATH_TO_DIR}/muri_metabarcoding/scripts
echo Config finished! $(date +"%T")
