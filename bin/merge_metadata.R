## Merging large sample metadata sheet and Miseq sample sheet
## 08/04/2023
## written by Ally Im

miseq_sheet <- read.csv("~/Desktop/dloop_sandbox/muri_metabarcoding/metadata/SampleSheetUsed.csv",header = TRUE, sep=",",skip = 19) #input sample sheet
sample_metadata <- read.csv("~/Desktop/dloop_sandbox/muri_metabarcoding/metadata/Hake_2019_metadata.csv",row.names = NULL) #input sample metadata file

miseq_sheet$sampleID <- paste0(sapply(strsplit(miseq_sheet$Sample_ID, "-"),`[`, 2),"-",sapply(strsplit(miseq_sheet$Sample_ID, "-"),`[`, 3))
print(miseq_sheet$sampleID)

metadata<- merge(x = miseq_sheet, y = sample_metadata, by="sampleID", all.x = TRUE)





