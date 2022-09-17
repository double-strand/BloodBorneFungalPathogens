## ---------------------------
## Script name: SILVA28SFungus.R
## Purpose of script: get Fungus 28S sequence data from latest SILVA resource   
## restriction site coordinates and flanking sequences.   
##
## Date Created: 2022-09-17
##
## Email: manuelxduval@double-strand.com
## ---------------------------
## Notes:
## 
## ---------------------------
# Dependencies.
sapply(c("Biostrings", "data.table", "dplyr", "ggplot2", 
         "kableExtra", "msa", "seqinr", "rvest","taxize","tidyverse"), 
       library, character.only = TRUE)
# setting the list of fungal sepcies to detect  
TargetSpecies <- c("Mucorales", "Aspergillus", "Fusarium", "Alternaria", 
                   "Scedosporium", "Aureobasidium", "Candida", "Cryptococcus", 
                   "Malassezia", "Pneumocystis", "Histoplasma", "Blastomyces", 
                   "Coccidioides", "Trichophyton", "Microsporum")

# getting SILVA database latest version number  
silva_latest_version <- str_trim(read_file("https://www.arb-silva.de/fileadmin/silva_databases/current/VERSION.txt"))

# reading the SILVA LSU 28S sequence data   
silva28S <- readDNAStringSet("https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz")

# selecting fungus entries
silvaFungi28S <- silva28S[grep("Fungi", names(silva28S), value = T)]

# setting SILVA species into a dataframe object
SILVA_28S_FungalSpecies <- as.data.frame(names(silvaFungi28S)) %>% tidyr::separate(`names(silvaFungi28S)`, c("SILVA_ID", "taxonomy"), sep = "\\s")



