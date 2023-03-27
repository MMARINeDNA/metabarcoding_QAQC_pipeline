#!/bin/sh

MFU_F="GTCGGTAAAACTCGTGCCAGC"
MFU_R="CATAGTGGGGTATCTAATCCCAGTTTG"
DL_F="TCACCCAAAGCTGRARTTCTA"
DL_R="GCGGGTTGCTGGTTTCACG"

#cd raw_fastqs
CUTADAPT=$(which cutadapt)
for i in *
do
FILE_PRIM=$(echo ${i} | cut -d - -f 1)
FILE_NAME=$(echo ${i} | cut -d . -f 1)
if [[ ${FILE_PRIM} == "MFU" ]]; then
if echo ${i} | grep -q "R1" ; then
echo MFU R1 Detected
${CUTADAPT} -a ${MFU_F} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
${i}> trimmed_${i} 2> ../cutadapt_reports/${FILE_NAME}_trim_report.txt
else
echo MFU R2 Detected
${CUTADAPT} -a ${MFU_R} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
${i}> trimmed_${i} 2> ../cutadapt_reports/${FILE_NAME}_trim_report.txt
fi
elif [[ ${FILE_PRIM} == "DL" ]]; then
if echo ${i} | grep -q "R1" ; then
echo DL R1 Detected
${CUTADAPT} -a ${DL_F} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
${i}> trimmed_${i} 2> ../cutadapt_reports/${FILE_NAME}_trim_report.txt
else
echo DL R2 Detected
${CUTADAPT} -a ${DL_R} \
    --discard-untrimmed \
    -j 0 \
    -q 30 \
${i}> trimmed_${i} 2> ../cutadapt_reports/${FILE_NAME}_trim_report.txt
fi
fi
done

cd ../cutadapt_reports
for i in *
do 
FILE_NAME=$(echo ${i} | cut -d . -f 1)
echo ${FILE_NAME} >> overall_report.txt 
grep Quality-trimmed $i >> overall_report.txt 
grep "Reads with adapters" $i >> overall_report.txt 
grep "Reads written (passing filters):" $i >> overall_report.txt 
done