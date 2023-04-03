## dada2 QAQC of MURI primer test sequences
## 11/28/2022
## written by Amy Van Cise using dada2

## set up working environment -------------------------------------------------------
library(dada2)
library(tidyverse)
library(seqinr)

args <- commandArgs(trailingOnly=TRUE)

fecal.seqs.file <- args[1]
  #" ~/Desktop/muri_sandbox/fastqs/"
taxref <- args[2]
  # " ~/Desktop/muri_sandbox/example_data_structure/metadata/cetacean_dloop_taxonomy.fasta"

### read primer test metadata ------------------------------------------------------

primer.data <- read.csv(args[3])
# ~/Desktop/muri_sandbox/example_data_structure/metadata/MURI301.csv
metadata_location <-"~/Desktop/muri_sandbox/example_data_structure/metadata/"

primer.data.pruned <- primer.data %>% 
  group_by(locus_shorthand) %>% 
  arrange(desc(primer_length)) %>% 
  slice_head()

#check if samples for i'th primer is present
for (i in 1:nrow(primer.data.pruned)){
  check <- grep(primer.data.pruned$locus_shorthand[i], sort(list.files(fecal.seqs.file, pattern="_R1_001.fastq", full.names = TRUE)), value = TRUE) 
  if(length(check) > 0){
  
      print(check)
### read fastq files in working directory -------------------------------------------
fnFs <- grep(primer.data.pruned$locus_shorthand[i], sort(list.files(fecal.seqs.file, pattern="_R1_001.fastq", full.names = TRUE)), value = TRUE)
fnRs <- grep(primer.data.pruned$locus_shorthand[i], sort(list.files(fecal.seqs.file, pattern="_R2_001.fastq", full.names = TRUE)), value = TRUE)
sample.names1 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names <- paste(sample.names1, sample.names2, sep = "_")
taxref <- grep(primer.data.pruned$locus_shorthand[i],list.files(path = metadata_location),value=TRUE)
tax_location <- paste0(metadata_location,taxref)
print("Running with Tax Database:")
print(tax_location)
### vizualize read quality profiles
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])


### Name filtered files in filtered/subdirectory ----------------------------------
filtFs <- file.path(fecal.seqs.file, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(fecal.seqs.file, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


### Filter and Trim ---------------------------------------------------------------

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                      trimLeft = 0, 
                     truncLen=0,
                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE)

# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
#                     trimLeft = primer.data.pruned$primer_length[i], 
#                      truncLen=c(ifelse(primer.data.pruned$tapestation_amplicon_length[i] > 300, 230, primer.data.pruned$tapestation_amplicon_length[i] - 50),
#                                 ifelse(primer.data.pruned$tapestation_amplicon_length[i] > 300, 230, primer.data.pruned$tapestation_amplicon_length[i] - 50)),
#                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                      compress=TRUE, multithread=TRUE)

### Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

### Learn Error Rates
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out #dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out

#plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)

### Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

### Merge Paired Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

### Construct sequence table
seqtab <- makeSequenceTable(mergers)

### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

freq.nochim <- sum(seqtab.nochim)/sum(seqtab)

### Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

### Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, tax_location, tryRC = TRUE, verbose = TRUE, multithread = TRUE, taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

### Save data
save(seqtab.nochim, freq.nochim, track, taxa, file = paste0("MURI_primer_test_mastertax_dada2_QAQC_output", primer.data.pruned$locus_shorthand[i], ".Rdata", sep = ""))
}
}


# use hashing for dada2!!!!
