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
TargetSpecies <- read_file("assets/targetSpecies.txt", header = T, skip = 2)

# getting SILVA database latest version number  
silva_latest_version <- str_trim(read_file("https://www.arb-silva.de/fileadmin/silva_databases/current/VERSION.txt"))

# SILVA URL
silva_url <- "https://www.arb-silva.de/fileadmin/silva_databases/"

# reading the SILVA LSU 28S sequence data   
silva28S <- readDNAStringSet(
  paste0(
    silva_url, 
    "release_", 
    silva_latest_version, 
    "/Exports//SILVA_", 
    silva_latest_version, 
    "_LSURef_NR99_tax_silva.fasta.gz"
    )
  )

# selecting fungus entries
silvaFungi28S <- silva28S[grep("Fungi", names(silva28S), value = T)]

# setting SILVA species into a dataframe object
SILVA_28S_FungalSpecies <- as.data.frame(names(silvaFungi28S)) %>%
  tidyr::separate(`names(silvaFungi28S)`, c("SILVA_ID", "taxonomy"), sep = "Eukaryota;")

# subsetting the SILVA_28S_Fungal Species for the target species 
SILVA_28S_TargetSpecies <- SILVA_28S_FungalSpecies[unlist(sapply(TargetSpecies, function(s) grep(s, SILVA_28S_FungalSpecies$taxonomy))),] %>% 
  tidyr::separate(., col = taxonomy, into = c("Clade", "Supergroup", "Group", "Kingdom", "Phylum", "Subphylum", "Class", "Order", "Family", "Genus", "Species", "Isolate"), ";", convert = T)








