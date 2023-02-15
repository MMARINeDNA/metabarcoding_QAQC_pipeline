## dada2 taxonomic classification of MURI primer test sequences
## 1/10/2023
## written by Amy Van Cise and Megan Shaffer using dada2

## set up working environment -------------------------------------------------------

library(dada2)
library(tidyverse)
library(seqinr)

taxref <- "H:/My Drive/01 MMARINeDNA MURI/01 Module 3/02 Bioinformatic analysis/01.5 Reference Taxonomy/02 reference fasta/MIDORI2_UNIQ_NUC_GB253_CO1_DADA2.fasta.gz"
loci.to.keep <- c("LRY")
taxonomy_name <- "MIDORI_CO1"

## Load QAQC'ed data -----------------------------------------------------------------

output.files <- list.files("H:/My Drive/01 MMARINeDNA MURI/01 Module 3/02 Bioinformatic analysis/02 pipeline output Rdata/", 
                           full.names = TRUE, pattern = "*.Rdata") 

# filter for the primers we want to use on this particular taxonomy
output.files <- output.files[which(grepl(paste(loci.to.keep, collapse = "|"),output.files))]


output.list <- list()

for (i in 1:length(output.files)){
  load(output.files[[i]],  temp_env <- new.env())
  output.list[[i]] <- as.list(temp_env)
}
rm(temp_env)


### Assign Taxonomy

for (i in 1:length(output.files)){

output.list[[i]]$MIDORI_C16_C18 <- assignTaxonomy(output.list[[i]]$seqtab.nochim, taxref, tryRC = TRUE)

temp <- output.list[[i]] 
save(temp, file = paste0("H:/My Drive/01 MMARINeDNA MURI/01 Module 3/02 Bioinformatic analysis/02 pipeline output Rdata/MURI_", loci.to.keep[i], "_", taxonomy_name, ".Rdata"))
rm(temp)

}

save(output.list, file = paste0("H:/My Drive/01 MMARINeDNA MURI/01 Module 3/02 Bioinformatic analysis/02 pipeline output Rdata/MURI_", paste(loci.to.keep, collapse = "_"), "_", taxonomy_name, ".Rdata"))
