# metabarcoding_QAQC_pipeline

This pipeline performs metabarcoding analysis using Cutadapt and DADA2.  

## Step 1: Pull the Repository
--- 
Before running this metabarcoding pipeline, it's important to make sure you have the most updated version of scripts and (more importantly) found ASV databases.  To do this, pull the most recent version of the github repo:

```
# run from the github repository
git pull
```
## Step 2: Create File System/Update Files
---

Below is the file system that this pipeline assumes:

![Alt text](./data/pictures/file_structure.png?raw=true "Title")

To create this filesystem on your computer, run the following line from the repo directory:

```
sh config.sh {path where to create files}
```
For example, if I wanted to create these files on my desktop, I would run:

```
sh config.sh ~/Desktop/
```
If you already have the filesystem on your system, the following command will copy over the updated files from the github repo to the filesytem

```
cp ./bin/* {path where to create files}/scripts
cp -r ./data/* {path where to create files}/metadata
```

## Step 3: Move Raw Fastq's into raw_fastq
---
Copy your raw fastq's into the raw_fastq directory.  A command that would do this would look like:

```
cp /path/to/fastqs/* /path/to/raw_fastqs
```

## Step 4: Run metabarcoding_wrapper.sh
The metabarcoding wrapper takes 2 inputs:
* path to your file system
* run name

To run this script, use the following command:
```
bash metabarcoding_wrapper.sh {pathway to files} {run name}
```

## Step 5: Push the updated ASV databases
To keep the ASV databases updated (and to make )

