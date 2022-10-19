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

primer4Footer <- "PRIMER_TASK=generic\nPRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_INTERNAL_OLIGO=1\nPRIMER_PICK_RIGHT_PRIMER=1\nPRIMER_OPT_SIZE=18\nPRIMER_MIN_SIZE=16\nPRIMER_MAX_SIZE=25\nPRIMER_MAX_NS_ACCEPTED=2\nPRIMER_PRODUCT_SIZE_RANGE=70-150\nP3_FILE_FLAG=1\nPRIMER_EXPLAIN_FLAG=1\n="

# setting the list of fungal sepcies to detect  
TargetSpecies <- read_table("assets/targetSpecies.txt", col_names = T, skip = 2)

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
# getting 28S sequence length distribution
SILVA28S_seqLength.stat <- summary(sapply(1:length(silvaFungi28S), function(i) length(silvaFungi28S[[i]])))
silvaFungi28S_seqLength <- sapply(1:length(silvaFungi28S), function(i) length(silvaFungi28S[[i]]))
#filtering out sequence entries that are either too long or too short
silvaFungi28S <- silvaFungi28S[which(silvaFungi28S_seqLength > SILVA28S_seqLength.stat[2] & silvaFungi28S_seqLength < SILVA28S_seqLength.stat[5])]
  

###############################################################################
# setting SILVA species into a dataframe object
SILVA_28S_FungalSpecies <- as.data.frame(names(silvaFungi28S)) %>%
  tidyr::separate(`names(silvaFungi28S)`, c("SILVA_ID", "taxonomy"), sep = "Eukaryota;")


# subsetting the SILVA_28S_Fungal Species for the target species 
SILVA_28S_TargetSpecies <- SILVA_28S_FungalSpecies[unlist(sapply(TargetSpecies$TaxonomicLabel, function(s) grep(s, SILVA_28S_FungalSpecies$taxonomy))),] %>% 
  tidyr::separate(., col = taxonomy, into = c("Clade", "Supergroup", "Group", "Kingdom", "Phylum", "Subphylum", "Class", "Order", "Family", "Genus", "Species", "Isolate"), ";", convert = T)

# setting a Biostrings object with the sequences of the target species
SILVA_28S_TargetSpecies.seq <- silvaFungi28S[names(silvaFungi28S) %in% unlist(lapply(TargetSpecies$TaxonomicLabel, function(s) grep(s, names(silvaFungi28S), value = T)))]

# selecting full length 28S rRNA seq
SILVA28S_seqLength.stat <- summary(sapply(1:length(SILVA_28S_TargetSpecies.seq), function(i) length(SILVA_28S_TargetSpecies.seq[[i]])))
SILVA28S_seqLength <- sapply(1:length(SILVA_28S_TargetSpecies.seq), function(i) length(SILVA_28S_TargetSpecies.seq[[i]]))

SILVA_28S_TargetSpecies.seq <- SILVA_28S_TargetSpecies.seq[which(SILVA28S_seqLength > SILVA28S_seqLength.stat[2] & SILVA28S_seqLength < SILVA28S_seqLength.stat[5])]

SILVA_28S_TargetSpecies.msa <- msa(SILVA_28S_TargetSpecies.seq[1:4], method = "ClustalW", gapOpening = 30, gapExtension = 5, maxiters = 20)
SILVA_28S_TargetSpecies.cons <- msaConsensusSequence(SILVA_28S_TargetSpecies.msa)
SILVA_28S_TargetSpecies.cons <- msaConsensusSequence(SILVA_28S_TargetSpecies.msa, thresh = c(80, 20), type="upperlower")
SILVA_28S_TargetSpecies.cons <- msaConsensusSequence(SILVA_28S_TargetSpecies.msa, thresh = c(100, 0), type="upperlower")

SILVA_28S_TargetSpecies.cons <- msaConsensusSequence(SILVA_28S_TargetSpecies.msa, thresh = c(90, 10), type="upperlower")

SILVA_28S_TargetSpecies.cons <- gsub("-|\\?", "N", SILVA_28S_TargetSpecies.cons)

SILVA_28S_TargetSpecies.cons <- DNAStringSet(SILVA_28S_TargetSpecies.cons)
names(SILVA_28S_TargetSpecies.cons) <- "SILVA_28S_TargetSpecies.cons"
writeXStringSet(SILVA_28S_TargetSpecies.cons, "SILVA_28S_TargetSpecies.cons")
SILVA_28S_TargetSpecies.cons <- read_templates("SILVA_28S_TargetSpecies.cons")

settings.xml <- system.file("extdata", "settings", "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
settings <- read_settings(settings.xml)
tmp <- design_primers(SILVA_28S_TargetSpecies.cons, mode.directionality = "fw", settings = design.settings)


SILVA_28S <- paste0("SEQUENCE_TEMPLATE=", SILVA_28S_TargetSpecies.cons)

primer4Header <- "SEQUENCE_ID=SILVA_138.1_LSURef_NR99"



fasta.file <- system.file("extdata", "IMGT_data", "templates",  "Homo_sapiens_IGH_functional_exon.fasta", package = "openPrimeR")
# Load the template sequences from 'fasta.file'
seq.df.simple <- read_templates(fasta.file)

################################################################################
####Mucormycetes
silvaChaetothyriales28S <- silvaFungi28S[grep("Chaetothyriales", names(silvaFungi28S), value = T)]

#shorten the name
names(silvaChaetothyriales28S) <- gsub("Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;Fungi;Dikarya;Ascomycota;Pezizomycotina", "", names(silvaChaetothyriales28S))

silvaChaetothyriales28S.msa <- msa(silvaChaetothyriales28S, method = "ClustalW", gapOpening = 50)
 
#silvaChaetothyriales28S.cons <- msaConsensusSequence(silvaChaetothyriales28S.msa, type = "upperlower")

silvaChaetothyriales28S.ali <- msaConvert(silvaChaetothyriales28S.msa, type = "seqinr::alignment")

silvaChaetothyriales28_dist <- seqinr::dist.alignment(silvaChaetothyriales28S.ali, matrix = "identity")

silvaChaetothyriales28.nj <- ape::nj(silvaChaetothyriales28_dist)

ape::plot.phylo(silvaChaetothyriales28.nj, main="Phylogenetic Tree", use.edge.length = F)

##Saccharomycetes
silvaSaccharomycetes28S <- silvaFungi28S[grep("Saccharomycetes", names(silvaFungi28S), value = T)]


silvaSaccharomycetes28S.msa <- msa(silvaSaccharomycetes28S, method = "ClustalW")

####
#seq1 <- unlist(strsplit(as.character(silvaSaccharomycetes28S[[1]]), ""))
TwentyMers <- kmer::kcount(seq1[1:30], k=20, residues = "DNA")
seq1 <- as.character(silvaSaccharomycetes28S[[1]])

# getting all 20-mers from the 1st sequence
Saccharomycetes28S_17meres <- sapply(0:2530, function(i) substring(as.character(silvaSaccharomycetes28S[[1]]), 1+i, 17+i))

#common20meres <- lapply(1:2530, function(i) grep(Saccharomycetes28S_20meres[i], as.character(silvaSaccharomycetes28S[[2]])))
#common20meres <- unlist(lapply(seq_along(common20meres), function(i) { rep(i, sum(common20meres[[i]] > 0)) }))

common20meres.list <- lapply(2:length(silvaSaccharomycetes28S), function(j) lapply(1:2530, function(i) grep(Saccharomycetes28S_20meres[i], as.character(silvaSaccharomycetes28S[[j]]))))

common20meres <- lapply(1:650, function(n) unlist(lapply(seq_along(common20meres.list[[n]]), function(i) { rep(i, sum(common20meres.list[[n]][[i]] > 0)) })))

# getting the overlapping cytometric quantities
common20meres.vec <- Reduce(intersect, sapply(1:3, function(i) common20meres[[i]]))




