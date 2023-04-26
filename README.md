# metabarcoding_QAQC_pipeline

This pipeline performs metabarcoding analysis using Cutadapt and DADA2.  Below is the file system that this pipeline assumes:

![Alt text](./data/pictures/file_structure.png?raw=true "Title")

To create this filesystem on your computer, run the following line from the repo directory:

```
sh config.sh {path where to create files}
```
For example, if I wanted to create these files on my desktop, I would run:

```
sh config.sh ~/Desktop/
```

