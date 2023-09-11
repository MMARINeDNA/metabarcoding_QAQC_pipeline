# MURI Metabarcoding Pipeline

This pipeline performs metabarcoding analysis using Cutadapt and DADA2. For more information on usage and dependencies, visit the wiki: https://github.com/MMARINeDNA/metabarcoding_QAQC_pipeline/wiki

For test data, visit the Zenodo: https://zenodo.org/record/8303276

## Overview of Metabarcoding Pipeline:

<p align="center">
<img src="https://github.com/MMARINeDNA/metabarcoding_QAQC_pipeline/blob/main/metadata/pictures/flowchart.png" alt="photo of metadata pipeline" width="400" class="center"/>
</p>

### **Step 1: Input Files** 
Take fastq files from Illumina sequencer and prepare data for pipeline 
### **Step 2: Primer Trimming** 
Use Cutadapt to trim primer sequences off of paired end sequences
### **Step 3: QC, Filtering, and Taxonomic Assignment**
Use DADA2 to perform quality trimming, filtering, paired-read merging, and more data manipulation.  Assign taxonomy using DADA2's naive bayes classifier.
### **Step 4: Preliminary Analysis**
Perform preliminary analysis with pyloseq and quarto.
### **Step 5: Output Files**
Review ouput files.



