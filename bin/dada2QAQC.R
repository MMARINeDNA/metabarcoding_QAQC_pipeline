
## dada2 QAQC of MURI primer test sequences
## 11/28/2022
## written by Amy Van Cise using dada2

## set up working environment -------------------------------------------------------
library(dada2)
library(tidyverse)
library(seqinr)
library(ShortRead)
library(digest)


### read in input files and variables ------------------------------------------------------
#args <- commandArgs(trailingOnly=TRUE)
fastq_location <- args[1]
output_location <- args[2]
metadata_location <- args[3]
run_name <- args[4]
primer.data <- read.csv(paste0(metadata_location,"primer_data.csv"))


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


### find taxonomy databases in metadata directory -------------------------------------------
    taxref <- grep(primer.data$locus_shorthand[i],list.files(path = metadata_location),value=TRUE)
    find_asv <- grep(primer.data$locus_shorthand[i],list.files(path = paste0(metadata_location,"known_hashes/")),value=TRUE)
    tax_location <- paste0(metadata_location,taxref)
    identified_hashes <- paste0(metadata_location,"known_hashes/",find_asv)
    print("Running with Tax Database:")
    print(tax_location)
    print("Running with ASV Database:")
    print(identified_hashes)


### Name filtered files in filtered/subdirectory ----------------------------------
    filtFs <- file.path(fastq_location, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(fastq_location, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names
    
### Filter out Empty Samples ----------------------------------
    file.empty <- function(filenames) file.info(filenames)$size == 20
    empty_files <- file.empty(fnFs) | file.empty(fnRs)
    fnFs <- fnFs[!empty_files]
    fnRs <- fnRs[!empty_files]
    filtFs <- filtFs[!empty_files]
    filtRs <- filtRs[!empty_files]
    sample.names <- sample.names[!empty_files]
    
### Find quality trimming length ----------------------------------
     print(paste0("Calculating quality trimming length...", Sys.time()))
     n <- 500000
     trimsF <- c()
     for(f in fnFs[!is.na(fnFs)]) {
       srqa <- qa(f, n=n)
       df <- srqa[["perCycle"]]$quality
       means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle) #calculate mean qual at each cycle
       where_to_cut <- min(which(means<30))-1 #trim first time mean qual dips below 30
       trimsF <- append(where_to_cut, trimsF) 
     }
     where_trim_all_Fs <- median(trimsF) #get average of all trims - use this as Trunclen forwards
     
     trimsR <- c()
     for(r in fnRs[!is.na(fnRs)]) {
       srqa <- qa(r, n=n)
       df <- srqa[["perCycle"]]$quality
       # Calculate summary statistics at each position
       means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)
       where_to_cut <- min(which(means<30))-1
       trimsR <- append(where_to_cut, trimsR)
     }
     where_trim_all_Rs <- median(trimsR)
     #try sliding window rule instead of hard cutoff
     both_length <- where_trim_all_Fs + where_trim_all_Rs
     min_length <- primer.data$tapestation_amplicon_length_F[i] + 25
     if(both_length < min_length){
       stop("Trim qual too high- not enough overlap. Choose a new Q score.") #if you trim too much, can't overlap
     }
     if(primer.data$locus_shorthand[i] == "DL"){
      if(where_trim_all_Fs > 280 || where_trim_all_Rs > 280){
        where_trim_all_Fs <- 280
        where_trim_all_Rs <- 280
      }
     }
    print(paste0("Finished calculating quality trimming length...", Sys.time()))
    
### Filter and Trim ---------------------------------------------------------------
    print(paste0("Starting filter and trim...", Sys.time()))
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                         trimRight = c(primer.data$primer_length_r[i],primer.data$primer_length_f[i]),
                         truncLen = c(150,150),
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                           compress=TRUE, multithread=TRUE, matchIDs=TRUE)
    print(paste0("Finished filter and trim.  ", Sys.time()))
    
### Dereplicate ---------------------------------------------------------------
    exists <- file.exists(filtFs) & file.exists(filtRs)
    filtFs <- filtFs[exists]
    filtRs <- filtRs[exists]
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    derepRs <- derepFastq(filtRs, verbose=TRUE)

    # Name the derep-class objects by the sample names
    sample.names <- sample.names[exists]
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names

### Learn Error Rates ---------------------------------------------------------------
    dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
    errF <- dadaFs.lrn[[1]]$err_out
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

### Filter by Size ---------------------------------------------------------------
    indexes.to.keep <- which(nchar(colnames(seqtab.nochim)) < primer.data$amplicon_length[i])
    cleaned.seqtab.nochim <- seqtab.nochim[,indexes.to.keep]
    
### Track reads through pipeline ---------------------------------------------------------------
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out[exists,], sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim),rowSums(cleaned.seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","length_filter")
    rownames(track) <- sample.names
    head(track)
    
### Create Hashing  ---------------------------------------------------------------
    
    # define output files
    conv_file <- file.path(output_location,"csv_output/",paste0(run_name,"_",primer.data$locus_shorthand[i],"_hash_key.csv"))
    conv_file.fasta <- file.path(output_location,"csv_output/",paste0(run_name,"_",primer.data$locus_shorthand[i],"_hash_key.fasta"))
    ASV_file <-  file.path(output_location,"csv_output/",paste0(run_name,"_",primer.data$locus_shorthand[i],"_ASV_table.csv"))
    taxonomy_file <- file.path(output_location,"csv_output/",paste0(run_name,"_",primer.data$locus_shorthand[i],"_taxonomy_output.csv"))
    
    # create ASV table and hash key 
    print(paste0("creating ASV table and hash key...", Sys.time()))
    seqtab.nochim.df <- as.data.frame(cleaned.seqtab.nochim)
    conv_table <- tibble( Hash = "", Sequence ="")
    Hashes <- map_chr (colnames(seqtab.nochim.df), ~ digest(.x, algo = "sha1", serialize = F, skip = "auto"))
    conv_table <- tibble (Hash = Hashes,
                          Sequence = colnames(seqtab.nochim.df))
    
    write_csv(conv_table, conv_file) # write hash key into a csv
    write.fasta(sequences = as.list(conv_table$Sequence), # write hash key into a fasta
                names     = as.list(conv_table$Hash),
                file.out = conv_file.fasta)
    sample.df <- tibble::rownames_to_column(seqtab.nochim.df,"Sample_name")
    sample.df <- data.frame(append(sample.df,c(Label=primer.data$locus_shorthand[i]), after = 1))
    current_asv <- bind_cols(sample.df %>%
                                    dplyr::select(Sample_name, Label),
                                  seqtab.nochim.df)
    current_asv <- current_asv %>%
      pivot_longer(cols = c(- Sample_name, - Label),
                   names_to = "Sequence",
                   values_to = "nReads") %>%
          filter(nReads > 0)
    
    current_asv <- merge(current_asv,conv_table, by="Sequence") %>%
      select(-Sequ)
  
  
    write_csv(current_asv, ASV_file) # write asv table into a csv

### Separate out ASVs that have already been classified -----------------------------------------
    identified_hashes <- read_csv(identified_hashes)
    Hash_key <- read_csv(conv_file)
    Hash <- Hash_key %>% 
      dplyr::select(Hash, Sequence) %>% 
      distinct()
    new_hashes_set <- anti_join(Hash, identified_hashes, by = c("Hash" = "Hash")) 
    already_identified_hashes_set <- inner_join(Hash,identified_hashes,by=c("Hash" = "Hash"))
    new.seqtab.nochim <- seqtab.nochim.df %>% 
      dplyr::select(new_hashes_set$Sequence) %>% 
                          as.matrix()
    
### Assign Taxonomy to new ASVs ---------------------------------------------------------------
    print(paste0("Starting Taxonomy Assignment at ", Sys.time()))
    new_taxa <- assignTaxonomy(new.seqtab.nochim,tax_location, tryRC = TRUE, verbose = TRUE, multithread = TRUE)
    print(paste0("Finished Taxonomy Assignment at ", Sys.time()," ."))
    
  
### Re-join Old and New ASVs --------------------------------------------------------------
    already_identified_hashes_set <- already_identified_hashes_set %>%
      select(c(Hash,Sequence,Kingdom,Phylum,Class,Order,Family,Genus,Species)) 
    
    ready_new_taxa <- new_taxa %>% 
      as.data.frame() %>%
      rownames_to_column("Sequence") %>%
      mutate(Hash = new_hashes_set$Hash) %>%
      relocate(Hash) 
    
    joined_old_new_taxa <- bind_rows(already_identified_hashes_set,ready_new_taxa)
    
### If new annotation to species level, add to previous effort -----------------------------------------
    to_add <- ready_new_taxa %>%
        filter(is.na(Species) == FALSE) %>%
        select(c(Hash,Kingdom,Phylum,Class,Order,Family,Genus,Species))
    updated_identified_hashes <- bind_rows(to_add,identified_hashes)
        
        
### Save data ---------------------------------------------------------------
    write.csv(joined_old_new_taxa,taxonomy_file) #write taxonomy csv
    write.csv(updated_identified_hashes,paste0(metadata_location,"known_hashes/",find_asv)) #write updated ASV database
    save(seqtab.nochim, freq.nochim, track, joined_old_new_taxa, file = paste0(output_location,"rdata_output/",run_name,"_dada2_output", primer.data$locus_shorthand[i], ".Rdata", sep = ""))
  }
}
  
