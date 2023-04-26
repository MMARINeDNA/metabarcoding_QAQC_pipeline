#!/bin/bash

#################### PRIMER VARIABLE ASSIGNMENT ####################
MFU_F="GCCGGTAAAACTCGTGCCAGC"
MFU_R="CATAGTGGGGTATCTAATCCCAGTTTG"
DL_F="TCACCCAAAGCTGRARTTCTA"
DL_R="GCGGGTTGCTGGTTTCACG"
MV1_F="CGTGCCAGCCACCGCG"
MV1_R="GGGTATCTAATCCYAGTTTG"
C16_F="GACGAGAAGACCCTAWTGAGCT"
C16_R="AAATTACGCTGTTATCCCT"
C18_F="CGCGGCGCTACATATTAGAC"
C16_R="GCACTTAACCGACCGTCGAC"


CUTADAPT=$(which cutadapt)
for i in *R1*
do
FILE_PRIM=$(echo ${i} | cut -d - -f 1) #grab primer name at beginning of file name
FILE_NAME=$(echo ${i} | cut -d _ -f 1,2,3) #grab everything before R1 (aka sample name)
R1=${i}
R2=$(echo ${i} | sed 's/R1/R2/g')
if [[ ${FILE_PRIM} == "MFU" ]]; then
echo MFU Detected
${CUTADAPT} -g ${MFU_F} \
     -G "${MFU_R}" \
     -o ../for_dada2/${R1} \
     -p ../for_dada2/${R2} \
    --discard-untrimmed \
    -j 0 \
"${R1}" "${R2}" 1> "../cutadapt_reports/${FILE_NAME}_trim_report.txt"
elif [[ ${FILE_PRIM} == "DL" ]]; then
echo DL Detected
${CUTADAPT} -g ${DL_F} \
     -G "${DL_R}" \
     -o ../for_dada2/${R1} \
     -p ../for_dada2/${R2} \
    --discard-untrimmed \
    -j 0 \
"${R1}" "${R2}" 1> "../cutadapt_reports/${FILE_NAME}_trim_report.txt"
elif [[ ${FILE_PRIM} == "MV1" ]]; then
echo MV1 Detected
${CUTADAPT} -g ${MV1_F} \
     -G "${MV1_R}" \
     -o ../for_dada2/${R1} \
     -p ../for_dada2/${R2} \
    --discard-untrimmed \
    -j 0 \
"${R1}" "${R2}" 1> "../cutadapt_reports/${FILE_NAME}_trim_report.txt"
elif [[ ${FILE_PRIM} == "C16" ]]; then
echo C16 Detected
${CUTADAPT} -g ${C16_F} \
     -G "${C16_R}" \
     -o ../for_dada2/${R1} \
     -p ../for_dada2/${R2} \
    --discard-untrimmed \
    -j 0 \
"${R1}" "${R2}" 1> "../cutadapt_reports/${FILE_NAME}_trim_report.txt"
elif [[ ${FILE_PRIM} == "C18" ]]; then
echo C18 Detected
${CUTADAPT} -g ${C18_F} \
     -G "${C18_R}" \
     -o ../for_dada2/${R1} \
     -p ../for_dada2/${R2} \
    --discard-untrimmed \
    -j 0 \
"${R1}" "${R2}" 1> "../cutadapt_reports/${FILE_NAME}_trim_report.txt"

fi
done