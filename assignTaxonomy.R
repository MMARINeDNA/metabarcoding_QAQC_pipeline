## dada2 taxonomic classification of MURI primer test sequences
## 1/10/2023
## written by Amy Van Cise and Megan Shaffer using dada2

## set up working environment -------------------------------------------------------

library(dada2)
library(tidyverse)
library(seqinr)

args <- commandArgs(trailingOnly=TRUE)

taxref <- args[1]
#taxref <- "~/Desktop/muri_sandbox/example_data_structure/metadata/cetacean_dloop_taxonomy.fasta"
loci.to.keep <- c(args[2])
taxonomy_name <- args[2]

## Load QAQC'ed data -----------------------------------------------------------------

output.files <- list.files("~/Desktop/muri_sandbox/example_data_structure/for_more_tax", 
                           full.names = TRUE, pattern = "*.Rdata") 

# filter for the primers we want to use on this particular taxonomy
output.files <- output.files[which(grepl(paste(loci.to.keep, collapse = "|"),output.files))]
print(output.files)


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
save(temp, file = paste0("~/Desktop/muri_sandbox/example_data_structure/final_data/", loci.to.keep[i], "_", taxonomy_name, ".Rdata"))
rm(temp)
}

save(output.list, file = paste0("~/Desktop/muri_sandbox/example_data_structure/final_data/", paste(loci.to.keep, collapse = "_"), "_", taxonomy_name, ".Rdata"))
