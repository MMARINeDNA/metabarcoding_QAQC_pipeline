## dada2 QAQC of MURI primer test sequences
## 11/28/2022
## written by Amy Van Cise using dada2

## set up working environment -------------------------------------------------------
library(dada2)
library(tidyverse)
library(seqinr)
library(ShortRead)

### read in fastq's ------------------------------------------------------
#args <- commandArgs(trailingOnly=TRUE)
#fastq_location <- args[1]
fastq_location <- "~/Desktop/muri_sandbox/example_data_structure/for_dada2"


### read primer metadata ------------------------------------------------------
primer.data <- read.csv("~/Desktop/muri_sandbox/example_data_structure/metadata/primer_data.csv")


### check if samples for i'th primer is present -------------------------------------------
for (i in 1:nrow(primer.data)){
  check <- grep(primer.data$locus_shorthand[i], sort(list.files(fastq_location, pattern="_R1_001.fastq", full.names = TRUE)), value = TRUE) 
  if(length(check) > 0){
    print(check)
    
    
### read fastq files in working directory -------------------------------------------
    fnFs <- grep(primer.data$locus_shorthand[i], sort(list.files(fastq_location, pattern="_R1_001.fastq", full.names = TRUE)), value = TRUE)
    fnRs <- grep(primer.data$locus_shorthand[i], sort(list.files(fastq_location, pattern="_R2_001.fastq", full.names = TRUE)), value = TRUE)
    sample.names1 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
    sample.names <- paste(sample.names1, sample.names2, sep = "_")


### find taxonomy database in metadata directory -------------------------------------------
    metadata_location <-"~/Desktop/muri_sandbox/example_data_structure/metadata/"
    taxref <- grep(primer.data$locus_shorthand[i],list.files(path = metadata_location),value=TRUE)
    tax_location <- paste0(metadata_location,taxref)
    print("Running with Tax Database:")
    print(tax_location)


### vizualize read quality profiles ----------------------------------
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])


### Name filtered files in filtered/subdirectory ----------------------------------
    filtFs <- file.path(fastq_location, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(fastq_location, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names
    
    
### Find quality trimming length ----------------------------------
    n <- 500000
    trimsF <- c()
    for(f in fnFs[!is.na(fnFs)]) {
      srqa <- qa(f, n=n)
      df <- srqa[["perCycle"]]$quality
      # Calculate summary statistics at each position
      means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)
      where_to_cut <- min(which(means<35))-1
      trimsF <- append(where_to_cut, trimsF)
    }
    where_trim_all_Fs <- mean(trimsF)
    
    trimsR <- c()
    for(r in fnRs[!is.na(fnRs)]) {
      srqa <- qa(r, n=n)
      df <- srqa[["perCycle"]]$quality
      # Calculate summary statistics at each position
      means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)
      where_to_cut <- min(which(means<35))-1
      trimsR <- append(where_to_cut, trimsR)
    }
    where_trim_all_Rs <- mean(trimsR)
    #try sliding window rule instead of hard cutoff
    #separate by F's and R's
    

### Filter and Trim ---------------------------------------------------------------
    print("starting filter and trim")
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                         trimRight = c(primer.data$primer_length_r[i],primer.data$primer_length_f[i]),
                         truncLen = c(where_trim_all_Fs,where_trim_all_Rs),
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                           compress=TRUE, multithread=FALSE)

    # out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
    #                     trimLeft = primer.data$primer_length[i], 
    #                      truncLen=c(ifelse(primer.data$tapestation_amplicon_length[i] > 300, 230, primer.data$tapestation_amplicon_length[i] - 50),
    #                                 ifelse(primer.data$tapestation_amplicon_length[i] > 300, 230, primer.data$tapestation_amplicon_length[i] - 50)),
    #                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
    #                      compress=TRUE, multithread=TRUE)

### Dereplicate ---------------------------------------------------------------
    exists <- file.exists(filtFs) & file.exists(filtRs)
    filtFs <- filtFs[exists]
    filtRs <- filtRs[exists]
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names

### Learn Error Rates ---------------------------------------------------------------
    dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
    errF <- dadaFs.lrn[[1]]$err_out #dadaFs.lrn[[1]]$err_out
    dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
    errR <- dadaRs.lrn[[1]]$err_out

### Sample Inference ---------------------------------------------------------------
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
    dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

### Merge Paired Reads ---------------------------------------------------------------
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

### Construct sequence table ---------------------------------------------------------------
    seqtab <- makeSequenceTable(mergers)

### Remove chimeras ---------------------------------------------------------------
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    freq.nochim <- sum(seqtab.nochim)/sum(seqtab)

### Track reads through pipeline ---------------------------------------------------------------
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names


### Assign Taxonomy ---------------------------------------------------------------
    print(paste0("Starting Taxonomy Assignment at ", Sys.time()))
    taxa <- assignTaxonomy(seqtab.nochim,tax_location, tryRC = TRUE, verbose = TRUE, multithread = TRUE)

### Save data ---------------------------------------------------------------
    save(seqtab.nochim, freq.nochim, track, taxref, file = paste0("MURI_primer_test_mastertax_dada2_QAQC_output", primer.data$locus_shorthand[i], ".Rdata", sep = ""))
  }
}
  

