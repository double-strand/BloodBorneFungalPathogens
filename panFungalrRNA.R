## ---------------------------
## Script name: panFungalrRNA.R
## Purpose of script: design primers and probes   
##
## Date Created: 2023-02-03
##
## Email: manuelxduval@double-strand.com
## ---------------------------
## Notes:
## 
## ---------------------------
# Dependencies.
sapply(c("Biostrings", "data.table", "dplyr", "ggplot2", 
         "kableExtra", "msa", "seqinr", "rvest","tidyverse"), 
       library, character.only = TRUE)

## 28S
silva28S <- readRDS("silva28S.rds")
# selecting fungus entries
silvaFungi28S <- silva28S[grep("Fungi", names(silva28S), value = T)]
# getting 28S sequence length distribution
SILVA28S_seqLength.stat <- summary(sapply(1:length(silvaFungi28S), function(i) length(silvaFungi28S[[i]])))
silvaFungi28S_seqLength <- sapply(1:length(silvaFungi28S), function(i) length(silvaFungi28S[[i]]))
#filtering out sequence entries that are either too long or too short
silvaFungi28S <- silvaFungi28S[which(silvaFungi28S_seqLength > SILVA28S_seqLength.stat[2] & silvaFungi28S_seqLength < SILVA28S_seqLength.stat[5])]
names(silvaFungi28S) <- gsub("Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;Fungi;", "", names(silvaFungi28S))

PanelA_B_silvaLSU_Species.df
## 18S
silva18S <- readRDS("silva18S.rds")
# selecting fungus entries
silvaFungi18S <- silva18S[grep("Fungi", names(silva18S), value = T)]
# getting 18S sequence length distribution
SILVA18S_seqLength.stat <- summary(sapply(1:length(silvaFungi18S), function(i) length(silvaFungi18S[[i]])))
silvaFungi18S_seqLength <- sapply(1:length(silvaFungi18S), function(i) length(silvaFungi18S[[i]]))
#filtering out sequence entries that are either too long or too short
silvaFungi18S <- silvaFungi18S[which(silvaFungi18S_seqLength > SILVA18S_seqLength.stat[2] & silvaFungi18S_seqLength < SILVA18S_seqLength.stat[5])]
names(silvaFungi18S) <- gsub("Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;Fungi;", "", names(silvaFungi18S))

PanelA_B_silvaSSU_Species.df
##########~PanelA~#######
###Mucormycetes 28S
silvaMucorales28S <- silvaFungi28S[grep("Mucorales", names(silvaFungi28S), value = T)]

silvaMucorales28S.msa <- msa(silvaMucorales28S, method = "ClustalOmega")
silvaMucorales28S.cons <- rprimer::consensusProfile(silvaMucorales28S.msa, ambiguityThreshold = 0.1)

#silvaMucorales28S.cons <- rprimer::consensusProfile(silvaMucorales28S.msa)

for (i in 1:5){
  tryCatch({
    print(i)
    silvaMucorales28S.reagents <- rprimer::designOligos(silvaMucorales28S.cons, maxDegeneracyPrimer = i, maxDegeneracyProbe = 2)
    silvaMucorales28S.assay <- designAssays(silvaMucorales28S.reagents)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

silvaMucorales28S.reagents <- rprimer::designOligos(silvaMucorales28S.cons, maxDegeneracyPrimer = 1, maxDegeneracyProbe = 2)
silvaMucorales28S.assay <- designAssays(silvaMucorales28S.reagents)
silvaMucorales28S.assay <- unique(silvaMucorales28S.assay[,c(1,2,3,6:9,20,21,24,25,26,27,38,39,44,45,46,47,60,61)])

###Mucormycetes 18S
silvaMucorales18S <- silvaFungi18S[grep("Mucorales", names(silvaFungi18S), value = T)]

silvaMucorales18S.df <- PanelA_B_silvaSSU_Species.df[PanelA_B_silvaSSU_Species.df$set == "Mucorales",]

MucoralesTargetFamily <- Cs(Mucor, Rhizopus, Rhizomucor, Cunninghamella, Lictheimia, Saksenaea, Syncephalastrum,Circinella, Apophysomyces, Cokeromyces, Zygomycetes, Rhizomucor, Apophysomyces, Saksenaea, Cokeromyces, Syncephalastrum) 

silvaMucorales18STargetSpecies.df <- silvaMucorales18S.df[silvaMucorales18S.df$Family %in% MucoralesTargetFamily,]

#silvaMucorales18S.msa <- msa(silvaMucorales18S, method = "ClustalOmega")
silvaMucorales18S.msa <- readRDS("silvaMucorales18S.msa.rds")
silvaMucorales18S.cons <- rprimer::consensusProfile(silvaMucorales18S.msa, ambiguityThreshold =0.025)

for (i in 1:5){
  tryCatch({
    print(i)
    silvaMucorales18S.reagents <- rprimer::designOligos(silvaMucorales18S.cons, maxDegeneracyPrimer = i, maxDegeneracyProbe = i)
    silvaMucorales18S.assay <-  rprimer::designAssays(silvaMucorales18S.reagents)
    if(exists("silvaMucorales18S.assay")) {
      break
    }
    print(i)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


##########~PanelA~#######
###Sordariomycetes 18S
Sordariomycetes18S <- silvaFungi18S[grep("Sordariomycetes", names(silvaFungi18S), value = T)]

Sordariomycetes_subset18S <- Sordariomycetes18S[names(Sordariomycetes18S) %in% grep("Fusarium|Purpureocillium|Acremonium|Scopulariopsis|Scedosporium|Lomentospora|Sporothrix", names(Sordariomycetes18S), value = T)]

Sordariomycetes_subset18S.msa <- msa(Sordariomycetes_subset18S, method = "ClustalOmega")
Sordariomycetes_subset18S.cons <- rprimer::consensusProfile(Sordariomycetes_subset18S.msa, ambiguityThreshold = 0)

Sordariomycetes18S.assay <- findAssay("Sordariomycetes", Sordariomycetes_subset18S.msa)
#Sordariomycetes18S.assay$fungi <- "Sordariomycetes"
#Sordariomycetes18S.assay <- Sordariomycetes18S.assay %>% dplyr::relocate(fungi)
  
##########~PanelA~#######

##########~Panel-B~#######
###Saccharomycetes 28S
Saccharomycetes28S <- silvaFungi28S[grep("Saccharomycetes", names(silvaFungi28S), value = T)]
Saccharomycete_subse28S <- Saccharomycetes28S[names(Saccharomycetes28S) %in% grep("Candida|Pichia|Meyerozyma|Debaryomyces|Cyberlindnera|Nakaseomyces", names(Saccharomycetes28S), value = T)]
Saccharomycete_subse28S.msa <- msa(Saccharomycete_subse28S, method = "ClustalOmega")
Saccharomycete28S.assay <- findAssay("Saccharomycetes", Saccharomycete_subse28S.msa)

Saccharomycete28S.assay2 <- as.data.frame(Saccharomycete28S.assay) %>% dplyr::group_by(iupacSequenceFwd) %>% dplyr::filter(gcContentPr == max(unlist(gcContentPr))) %>% distinct


######
findAssay <- function(fungi, msa) {
  for (k in seq(0,.2,.025)) {
    #print(k)
    fungi.cons <- rprimer::consensusProfile(msa, ambiguityThreshold =k)
    for (i in 1:5){
      tryCatch({
        #print(i)
       fungi.reagents <- rprimer::designOligos(fungi.cons, maxDegeneracyPrimer = i, maxDegeneracyProbe = i)
        fungi.assay <-  rprimer::designAssays(fungi.reagents)
        if(exists("fungi.assay")) {
          break
        }
      },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    if(exists("fungi.assay")) {
      print(k)
      print(i)
      return(cbind(as.data.frame(fungi), unique(fungi.assay[,c(1,2,3,6:9,20,21,24,25,26,27,38,39,44,45,46,47,60,61)])))
      break
    }
    }
}

TargetSpecies <- c("Mucorales", "Eurotiales", "Sordariomycetes", "Dothideomycetes", "Saccharomycetes", "Basidiomycota", "Taphrinomycotina", "Onygenales", "Chaetothyriales")

panFungalFunction <- function(su, fungi) {
  allowed_su <- c("18S", "28S")
  if (!su %in% allowed_su) {
    return("Invalid input rna, select 18S or 28S")
  }
  allowed_fungi <- c("Mucorales", "Eurotiales", "Sordariomycetes", "Dothideomycetes", "Saccharomycetes", "Basidiomycota", "Taphrinomycotina", "Onygenales", "Chaetothyriales")
  if (!fungi %in% allowed_fungi) {
    return("Invalid fungi, select 18S or 28S")
  }
  if (su == "18S") {
    targetedFungi <- silvaFungi18S[grep(fungi, names(silvaFungi18S), value = T)]
  } else {
    targetedFungi <- silvaFungi28S[grep(fungi, names(silvaFungi28S), value = T)]
  }
  targetedFungi.msa <- msa(targetedFungi, method = "ClustalOmega")
  
  #findAssay <- function() {
    for (k in seq(0,.2,.025)) {
      #print(k)
      targetedFungi.cons <- rprimer::consensusProfile(targetedFungi.msa, ambiguityThreshold =k)
      for (i in 1:5){
        tryCatch({
          #print(i)
          targetedFungi.reagents <- rprimer::designOligos(targetedFungi.cons, maxDegeneracyPrimer = i, maxDegeneracyProbe = i)
          targetedFungi.assay <-  rprimer::designAssays(targetedFungi.reagents)
          if(exists("silvaMucorales18S.assay")) {
            break
          }
        },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
      if(exists("targetedFungi.assay")) {
        print(k)
        print(i)
        return(cbind(as.data.frame(fungi), unique(targetedFungi.assay[,c(1,2,3,6:9,20,21,24,25,26,27,38,39,44,45,46,47,60,61)])))
        break
      }
    }
  }
  #return(targetedFungi.cons)
#}
#Sordariomycetes
tmp <- panFungalFunction("18S", "Sordariomycetes")
silvaMucorales18S.msa <- readRDS("silvaMucorales18S.msa.rds")
silvaMucorales18S.cons <- rprimer::consensusProfile(silvaMucorales18S.msa, ambiguityThreshold =0.025)


Assay28S <- lapply(TargetSpecies, function(s) panFungalFunction("28S", s))
Assay28S.df <- do.call(rbind, Assay28S)

mucorales.df <- Assay28S.df[Assay28S.df$fungi == "Mucorales",]

table(mucorales.df$start)
#1023 1024 1025 1026 
#351  351  234  117 

mucorales.df1 <- mucorales.df[mucorales.df$start == 1023,]


Assay28S.df <- rbind(Assay28S.df, Saccharomycete28S.assay)
Assay28S.df$rRNA <- "28S"
Assay28S.df <- Assay28S.df %>% dplyr::relocate(rRNA)

saveRDS(Assay28S.df, "Assay28S.rds")

Assay18S <- lapply(TargetSpecies, function(s) panFungalFunction("18S", s))
Assay18S.df <- do.call(rbind, Assay18S)

Assay18S.df <- rbind(Assay18S.df, Sordariomycetes18S.assay)
Assay18S.df <- rbind(Assay18S.df, Dothideomycete28S.assay)
Assay18S.df <- rbind(Assay18S.df, Basidiomycota18S.assay)
Assay18S.df$rRNA <- "18S"

Assay18S.df <- Assay18S.df %>% dplyr::relocate(rRNA)
saveRDS(Assay18S.df, "Assay18S.rds")

panelA_BAssays.df <- rbind(Assay28S.df, Assay18S.df)

write_csv2(Assay18S.df, "panFungalAssays.csv")

#Sordariomycetes
#SordariomycetesAssays28S <- lapply("Sordariomycetes", function(s) panFungalFunction("28S", s))
Sordariomycetes28S <- silvaFungi28S[grep("Sordariomycetes", names(silvaFungi28S), value = T)]
Sordariomycetes_subset28S <- Sordariomycetes28S[names(Sordariomycetes28S) %in% grep("Fusarium|Purpureocillium|Acremonium|Scopulariopsis|Scedosporium|Lomentospora|Sporothrix", names(Sordariomycetes28S), value = T)]
Sordariomycetes_subset28S.msa <- msa(Sordariomycetes_subset28S, method = "ClustalOmega")
Sordariomycetes28S.assay <- findAssay("Sordariomycetes", Sordariomycetes_subset28S.msa)

Sordariomycetes18S <- silvaFungi18S[grep("Sordariomycetes", names(silvaFungi18S), value = T)]
Sordariomycetes_subset18S <- Sordariomycetes18S[names(Sordariomycetes18S) %in% grep("Fusarium|Purpureocillium|Acremonium|Scopulariopsis|Scedosporium|Lomentospora|Sporothrix", names(Sordariomycetes18S), value = T)]
Sordariomycetes_subset18S.msa <- msa(Sordariomycetes_subset18S, method = "ClustalOmega")
Sordariomycetes18S.assay <- findAssay("Sordariomycetes", Sordariomycetes_subset18S.msa)

Basidiomycota18S <- silvaFungi18S[grep("Basidiomycota", names(silvaFungi18S), value = T)]
Basidiomycota_subset18S <- Basidiomycota18S[names(Basidiomycota18S) %in% grep("Cryptococcus|Trichosporon|Malassezia", names(Basidiomycota18S), value = T)]
Basidiomycota_subset18S.msa <- msa(Basidiomycota_subset18S, method = "ClustalOmega")
Basidiomycota18S.assay <- findAssay("Basidiomycota", Basidiomycota_subset18S.msa)

Basidiomycota28S <- silvaFungi28S[grep("Basidiomycota", names(silvaFungi28S), value = T)]
Basidiomycota_subset28S <- Basidiomycota28S[names(Basidiomycota28S) %in% grep("Cryptococcus|Trichosporon|Malassezia", names(Basidiomycota28S), value = T)]
Basidiomycota_subset28S.msa <- msa(Basidiomycota_subset28S, method = "ClustalOmega")
Basidiomycota28S.assay <- findAssay("Basidiomycota", Basidiomycota_subset28S.msa)
Basidiomycota28S.assay2 <- as.data.frame(Basidiomycota28S.assay) %>% dplyr::group_by(iupacSequenceFwd) %>% dplyr::filter(gcContentPr == max(unlist(gcContentPr))) %>% distinct
Basidiomycota28S.assay3 <- as.data.frame(Basidiomycota28S.assay2) %>% dplyr::group_by(iupacSequenceFwd) %>% dplyr::filter(length == min(length)) %>% distinct


df <- read.table("../panfungalAssays/panfungalAssaySelection/assaysv1.csv", sep = "\t", header = T)
df <- rbind(df, Sordariomycetes18S.assay, Dothideomycete18S.assay[1:50,], Saccharomycete28S.assay2, Basidiomycota28S.assay3[1:50,])
write_tsv(df, "../panfungalAssays/panfungalAssaySelection/assays.csv")

###
#Ascomycota
# 28S
Saccharomycetes28S <- silvaFungi28S[grep("Saccharomycetes", names(silvaFungi28S), value = T)]
Saccharomycetes_subset28S <- Saccharomycetes28S[names(Saccharomycetes28S) %in% grep("Candida|Pichia|Meyerozyma|Debaryomyces|Cyberlindnera|Nakaseomyces", names(Saccharomycetes28S), value = T)]
Saccharomycetes_subset28S.msa <- msa(Saccharomycetes_subset28S, method = "ClustalOmega")
Saccharomycetes28S.assay <- findAssay("Saccharomycetes", Saccharomycetes_subset28S.msa)
#~ assays found for Saccharomycetes28S

Saccharomycetes18S <- silvaFungi18S[grep("Saccharomycetes", names(silvaFungi18S), value = T)]
Saccharomycetes_subset18S <- Saccharomycetes18S[names(Saccharomycetes18S) %in% grep("Candida|Pichia|Meyerozyma|Debaryomyces|Cyberlindnera|Nakaseomyces", names(Saccharomycetes18S), value = T)]
Saccharomycetes_subset18S.msa <- msa(Saccharomycetes_subset18S, method = "ClustalOmega")
Saccharomycetes18S.assay <- findAssay("Saccharomycetes", Saccharomycetes_subset18S.msa)
#~ No assays found for Saccharomycetes18S

Eurotiales28S <- silvaFungi28S[grep("Eurotiales", names(silvaFungi28S), value = T)]
Eurotiales_subset28S <- Eurotiales28S[names(Eurotiales28S) %in% grep("Aspergillus|Paecilomyces|Penicillium|Talaromyces", names(Eurotiales28S), value = T)]
Eurotiales_subset28S.msa <- msa(Eurotiales_subset28S, method = "ClustalOmega")
Eurotiales28S.assay <- findAssay("Eurotiales", Eurotiales_subset28S.msa)
#~ assays found for Eurotiales28S

Eurotiales18S <- silvaFungi18S[grep("Eurotiales", names(silvaFungi18S), value = T)]
Eurotiales_subset18S <- Eurotiales18S[names(Eurotiales18S) %in% grep("Aspergillus|Paecilomyces|Penicillium|Talaromyces", names(Eurotiales18S), value = T)]
Eurotiales_subset18S.msa <- msa(Eurotiales_subset18S, method = "ClustalOmega")
Eurotiales18S.assay <- findAssay("Eurotiales", Eurotiales_subset18S.msa)
#~ No assays found for Eurotiales18S

Chaetothyriales28S <- silvaFungi28S[grep("Chaetothyriales", names(silvaFungi28S), value = T)]
#~ only one entry for Chaetothyriales 28S

Chaetothyriales18S <- silvaFungi18S[grep("Chaetothyriales", names(silvaFungi18S), value = T)]
Chaetothyriales_subset18S <- Chaetothyriales18S[names(Chaetothyriales18S) %in% grep("Fonsecaea|Phialophora|Exophiala", names(Chaetothyriales18S), value = T)]
Chaetothyriales_subset18S.msa <- msa(Chaetothyriales_subset18S, method = "ClustalOmega")
Chaetothyriales18S.assay <- findAssay("Chaetothyriales", Chaetothyriales_subset18S.msa)
#~ assays found for Chaetothyriales18S

Onygenales28S <- silvaFungi28S[grep("Onygenales", names(silvaFungi28S), value = T)]
Onygenales_subset28S <- Onygenales28S[names(Onygenales28S) %in% grep("Chrysosporium|Histoplasma|Blastomyces|Coccidioides|Paracoccidioides|Microsporum|Epidermophyton|Trichophyton", names(Onygenales28S), value = T)]
Onygenales_subset28S.msa <- msa(Onygenales_subset28S, method = "ClustalOmega")
Onygenales28S.assay <- findAssay("Onygenales", Onygenales_subset28S.msa)
#~ no assays found for Onygenales28S

Onygenales18S <- silvaFungi18S[grep("Onygenales", names(silvaFungi18S), value = T)]
Onygenales_subset18S <- Onygenales18S[names(Onygenales18S) %in% grep("Chrysosporium|Histoplasma|Blastomyces|Coccidioides|Paracoccidioides|Microsporum|Epidermophyton|Trichophyton", names(Onygenales18S), value = T)]
Onygenales_subset18S.msa <- msa(Onygenales_subset18S, method = "ClustalOmega")
Onygenales18S.assay <- findAssay("Onygenales", Onygenales_subset18S.msa)
#~ assays found for Onygenales18S

Sordariomycetes28S <- silvaFungi28S[grep("Sordariomycetes", names(silvaFungi28S), value = T)]
Sordariomycetes_subset28S <- Sordariomycetes28S[names(Sordariomycetes28S) %in% grep("Fusarium|Purpureocillium|Acremonium|Scopulariopsis|Scedosporium|Lomentospora|Sporothrix", names(Sordariomycetes28S), value = T)]
Sordariomycetes_subset28S.msa <- msa(Sordariomycetes_subset28S, method = "ClustalOmega")
Sordariomycetes28S.assay <- findAssay("Sordariomycetes", Sordariomycetes_subset28S.msa)
#~ no assays found for Sordariomycetes 18S

Sordariomycetes18S <- silvaFungi18S[grep("Sordariomycetes", names(silvaFungi18S), value = T)]
Sordariomycetes_subset18S <- Sordariomycetes18S[names(Sordariomycetes18S) %in% grep("Fusarium|Purpureocillium|Acremonium|Scopulariopsis|Scedosporium|Lomentospora|Sporothrix", names(Sordariomycetes18S), value = T)]
Sordariomycetes_subset18S.msa <- msa(Sordariomycetes_subset18S, method = "ClustalOmega")
Sordariomycetes18S.assay <- findAssay("Sordariomycetes", Sordariomycetes_subset18S.msa)
#~ assays found for Sordariomycetes 18S

###Dothideomycetes 18S
Dothideomycete28S <- silvaFungi28S[grep("Dothideomycete", names(silvaFungi28S), value = T)]
Dothideomycete_subset28S <- Dothideomycete28S[names(Dothideomycete28S) %in% grep("Pleosporales|Cladosporiales|Dothideales", names(Dothideomycete28S), value = T)]
Dothideomycete_subset28S.msa <- msa(Dothideomycete_subset28S, method = "ClustalOmega")
Dothideomycete28S.assay <- findAssay("Dothideomycete", Dothideomycete_subset28S.msa)
#~ No assays found for Dothideomycete 28S

Dothideomycetes18S <- silvaFungi18S[grep("Dothideomycetes", names(silvaFungi18S), value = T)]
Dothideomycetes_subset18S <- Dothideomycetes18S[names(Dothideomycetes18S) %in% grep("Alternaria|Curvularia|Bipolaris|Exserohilum|Cladosporium|Aureobasidium", names(Dothideomycetes18S), value = T)]
Dothideomycetes_subset18S.msa <- msa(Dothideomycetes_subset18S, method = "ClustalOmega")
Dothideomycetes18S.assay <- findAssay("Dothideomycetes", Dothideomycetes_subset18S.msa)
#~ assays found for Dothideomycete 18S

Taphrinomycotina28S <- silvaFungi28S[grep("Taphrinomycotina", names(silvaFungi28S), value = T)]
Taphrinomycotina_subset28S <- Taphrinomycotina28S[names(Taphrinomycotina28S) %in% grep("Pneumocystis", names(Taphrinomycotina28S), value = T)]
#~ one entry

Taphrinomycotina18S <- silvaFungi18S[grep("Taphrinomycotina", names(silvaFungi18S), value = T)]
Taphrinomycotina_subset18S <- Taphrinomycotina18S[names(Taphrinomycotina18S) %in% grep("Pneumocystis", names(Taphrinomycotina18S), value = T)]
Taphrinomycotina_subset18S.msa <- msa(Taphrinomycotina_subset18S, method = "ClustalOmega")
Taphrinomycotina18S.assay <- findAssay("Taphrinomycotina", Taphrinomycotina_subset18S.msa)
#~ assays found for Taphrinomycotina 18S

# Basidiomycota
###Basidiomycota 18S
Basidiomycota28S <- silvaFungi28S[grep("Basidiomycota", names(silvaFungi28S), value = T)]
Basidiomycota_subset28S <- Basidiomycota28S[names(Basidiomycota28S) %in% grep("Cryptococcus|Trichosporon|Malassezia", names(Basidiomycota28S), value = T)]
Basidiomycota_subset28S.msa <- msa(Basidiomycota_subset28S, method = "ClustalOmega")
Basidiomycota28S.assay <- findAssay("Basidiomycota", Basidiomycota_subset28S.msa)
#~ assays found for Basidiomycota 28S

Basidiomycota18S <- silvaFungi28S[grep("Basidiomycota", names(silvaFungi28S), value = T)]
Basidiomycota_subset18S <- Basidiomycota18S[names(Basidiomycota18S) %in% grep("Cryptococcus|Trichosporon|Malassezia", names(Basidiomycota18S), value = T)]
Basidiomycota_subset18S.msa <- msa(Basidiomycota_subset18S, method = "ClustalOmega")
Basidiomycota18S.assay <- findAssay("Basidiomycota", Basidiomycota_subset18S.msa)
#~ assays found for Basidiomycota 18S

# Mucorales
Mucorales28S <- silvaFungi28S[grep("Mucorales", names(silvaFungi28S), value = T)]
Mucorales_subset28S <- Mucorales28S[names(Mucorales28S) %in% grep("Mucor|Rhizopus|Rhizomucor|Cunninghamella|Lictheimia|Saksenaea|Syncephalastrum|Circinella|Apophysomyces|Cokeromyces", names(Mucorales28S), value = T)]
Mucorales_subset28S.msa <- msa(Mucorales_subset28S, method = "ClustalOmega")
Mucorales_subset28S.assay <- findAssay("Mucorales", Mucorales_subset28S.msa)
#~ assays found for Mucorales 28S

Mucorales18S <- silvaFungi18S[grep("Mucorales", names(silvaFungi18S), value = T)]
Mucorales_subset18S <- Mucorales18S[names(Mucorales18S) %in% grep("Mucor|Rhizopus|Rhizomucor|Cunninghamella|Lictheimia|Saksenaea|Syncephalastrum|Circinella|Apophysomyces|Cokeromyces", names(Mucorales18S), value = T)]
Mucorales_subset18S.msa <- msa(Mucorales_subset18S, method = "ClustalOmega")
Mucorales_subset18S.assay <- findAssay("Mucorales", Mucorales_subset18S.msa)
#~ assays found for Mucorales 18S


#Onygenales and Sordariomycetes
Onygenales_Sordariomycetes_18S <- c(Onygenales_subset18S, Sordariomycetes_subset18S)
Onygenales_Sordariomycetes_18S.msa <- msa(Onygenales_Sordariomycetes_18S, method = "ClustalOmega")
Onygenales_Sordariomycetes_18S.assay <- findAssay("Onygenales_Sordariomycetes", Onygenales_Sordariomycetes_18S.msa)
#~ No assays found for nygenales and Sordariomycetes

#Taphrinomycotina and Chaetothyriales
Taphrinomycotina_Chaetothyriales_28S <- c(Taphrinomycotina_subset28S, Chaetothyriales28S)
Taphrinomycotina_Chaetothyriales_28S.msa <- msa(Taphrinomycotina_Chaetothyriales_28S,  method = "ClustalOmega")
Taphrinomycotina_Chaetothyriales_28S.assay <- findAssay("Taphrinomycotina_Chaetothyriale", Taphrinomycotina_Chaetothyriales_28S.msa)
#~ No assays found for Taphrinomycotina and Chaetothyriales28S
Taphrinomycotina_Chaetothyriales_18S <- c(Taphrinomycotina_subset18S, Chaetothyriales_subset18S)
Taphrinomycotina_Chaetothyriales_18S.msa <- msa(Taphrinomycotina_Chaetothyriales_18S,  method = "ClustalOmega")
Taphrinomycotina_Chaetothyriales_18S.assay <- findAssay("Taphrinomycotina_Chaetothyriale", Taphrinomycotina_Chaetothyriales_18S.msa)
#~ assays found for Taphrinomycotina and Chaetothyriales 18S

# Chaetothyriales, Onygenales, Taphrinomycotina and Sordariomycetes
Taphri_Chaetot_Onygenales_Sordariomycetes_18S <- c(Taphrinomycotina_subset18S, Chaetothyriales_subset18S, Onygenales_subset18S, Sordariomycetes_subset18S)
Taphri_Chaetot_Onygenales_Sordariomycetes_18S.msa <- msa(Taphri_Chaetot_Onygenales_Sordariomycetes_18S,  method = "ClustalOmega")
Taphri_Chaetot_Onygenales_Sordariomycetes_18S.assay <- findAssay("Taphrinomycotina_Chaetothyriales_Onygenales_Sordariomycetes", Taphri_Chaetot_Onygenales_Sordariomycetes_18S.msa)
# No assay

# Chaetothyriales, Onygenales, Taphrinomycotina and Dothideomycetes
Taphri_Chaetot_Onygenales_Dothideomycetes_18S <- c(Taphrinomycotina_subset18S, Chaetothyriales_subset18S, Onygenales_subset18S, Dothideomycetes_subset18S)
Taphri_Chaetot_Onygenales_Dothideomycetes_18S.msa <- msa(Taphri_Chaetot_Onygenales_Dothideomycetes_18S,  method = "ClustalOmega")
Taphri_Chaetot_Onygenales_Dothideomycetes_18S.assay <- findAssay("Taphrinomycotina_Chaetothyriales_Onygenales_Dothideomycetes", Taphri_Chaetot_Onygenales_Dothideomycetes_18S.msa)
# No assay

# Chaetothyriales, Onygenales, Taphrinomycotina and Basidiomycota
Taphri_Chaetot_Onygenales_Basidiomycota_18S <- c(Taphrinomycotina_subset18S, Chaetothyriales_subset18S, Onygenales_subset18S, Basidiomycota_subset18S)
Taphri_Chaetot_Onygenales_Basidiomycota_18S.msa <- msa(Taphri_Chaetot_Onygenales_Basidiomycota_18S,  method = "ClustalOmega")
Taphri_Chaetot_Onygenales_Basidiomycota_18S.assay <- findAssay("Taphrinomycotina_Chaetothyriales_Onygenales_Basidiomycota", Taphri_Chaetot_Onygenales_Basidiomycota_18S.msa)
# No assay



########## assays ###############################
#Taphrinomycotina and Chaetothyriales and Onygenales
Taphrinomycotina_Chaetothyriales_Onygenales_18S <- c(Taphrinomycotina_subset18S, Chaetothyriales_subset18S, Onygenales_subset18S)
Taphrinomycotina_Chaetothyriales_Onygenales_18S.msa <- msa(Taphrinomycotina_Chaetothyriales_Onygenales_18S, method = "ClustalOmega")
Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay <- findAssay("Taphrinomycotina_Chaetothyriales_Onygenales",Taphrinomycotina_Chaetothyriales_Onygenales_18S.msa)
Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay <- tidyr::unnest(Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay)
Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay <- Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay[Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay$tmPr >60,]
Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay <- Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay[Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay$length <73,]
# assays

# Saccharomycetales and Basidiomycota
Saccharomycetales_Basidiomycota_28S <- c(Saccharomycete_subset28S, Basidiomycota_subset28S)
Saccharomycetales_Basidiomycota_28S.msa <- msa(Saccharomycetales_Basidiomycota_28S, method = "ClustalOmega")
Saccharomycetales_Basidiomycota_28S.assay <- findAssay("Saccharomycetales_Basidiomycota", Saccharomycetales_Basidiomycota_28S.msa)
Saccharomycetales_Basidiomycota_28S.assay <- tidyr::unnest(Saccharomycetales_Basidiomycota_28S.assay)
Saccharomycetales_Basidiomycota_28S.assay<- Saccharomycetales_Basidiomycota_28S.assay[Saccharomycetales_Basidiomycota_28S.assay$tmPr > 56,]
Saccharomycetales_Basidiomycota_28S.assay<- Saccharomycetales_Basidiomycota_28S.assay[Saccharomycetales_Basidiomycota_28S.assay$length < 117,]
# assays 

# Sordariomycetes and Dothideomycetes
Sordariomycetes_Dothideomycetes_18S <- c(Sordariomycetes_subset18S, Dothideomycetes_subset18S)
Sordariomycetes_Dothideomycetes_18S.msa <- msa(Sordariomycetes_Dothideomycetes_18S, method = "ClustalOmega")
Sordariomycetes_Dothideomycetes_18S.assay <- findAssay("Sordariomycetes_Dothideomycetes", Sordariomycetes_Dothideomycetes_18S.msa)
Sordariomycetes_Dothideomycetes_18S.assay <- tidyr::unnest(Sordariomycetes_Dothideomycetes_18S.assay)
Sordariomycetes_Dothideomycetes_18S.assay <- Sordariomycetes_Dothideomycetes_18S.assay[Sordariomycetes_Dothideomycetes_18S.assay$tmPr > 65,]
Sordariomycetes_Dothideomycetes_18S.assay <- Sordariomycetes_Dothideomycetes_18S.assay[Sordariomycetes_Dothideomycetes_18S.assay$length < 114,]
# assays

# Mucorales and Eurotiales
Mucorales_Eurotiales_28S <- c(Mucorales_subset28S, Eurotiales_subset28S)
Mucorales_Eurotiales_28S.msa <- msa(Mucorales_Eurotiales_28S, method = "ClustalOmega")
Mucorales_Eurotiales_28S.assay <- findAssay("Mucorales_Eurotiales", Mucorales_Eurotiales_28S.msa)
# assays

Mucorales_Eurotiales_18S <- c(Mucorales_subset18S, Eurotiales_subset18S)
Mucorales_Eurotiales_18S.msa <- msa(Mucorales_Eurotiales_18S, method = "ClustalOmega")
Mucorales_Eurotiales_18S.assay <- findAssay("Mucorales_Eurotiales", Mucorales_Eurotiales_18S.msa)
Mucorales_Eurotiales_18S.assay <- tidyr::unnest(Mucorales_Eurotiales_18S.assay)
Mucorales_Eurotiales_18S.assay <- Mucorales_Eurotiales_18S.assay[Mucorales_Eurotiales_18S.assay$tmPr > 63,]
Mucorales_Eurotiales_18S.assay <- Mucorales_Eurotiales_18S.assay[Mucorales_Eurotiales_18S.assay$length > 71,]
# assays 

#############
# Sordariomycetes_Dothideomycetes_Mucorales_Eurotiales_18S <- c(Mucorales_subset18S, Eurotiales_subset18S, Sordariomycetes_subset18S, Dothideomycetes_subset18S)
# Sordariomycetes_Dothideomycetes_Mucorales_Eurotiales_18S.msa <- msa(Sordariomycetes_Dothideomycetes_Mucorales_Eurotiales_18S, method = "ClustalOmega")
# Sordariomycetes_Dothideomycetes_Mucorales_Eurotiales_18S.assay <- findAssay("Sordariomycetes_Dothideomycetes_Mucorales_Eurotiales", Sordariomycetes_Dothideomycetes_Mucorales_Eurotiales_18S.msa)
# # no assays
# 
# Saccharomycete_Basidiomycota_Mucorales_Eurotiales_18S <- c(Mucorales_subset28S, Eurotiales_subset28S, Saccharomycetes_subset28S, Basidiomycota_subset28S)
# Saccharomycete_Basidiomycota_Mucorales_Eurotiales_18S.msa <- msa(Saccharomycete_Basidiomycota_Mucorales_Eurotiales_18S, method = "ClustalOmega")
# Saccharomycete_Basidiomycota_Mucorales_Eurotiales_18S.assay <- findAssay("Saccharomycete_Basidiomycota_Mucorales_Eurotiales", Saccharomycete_Basidiomycota_Mucorales_Eurotiales_18S.msa)
# # no assays
# 
# Sordariomycetes_Dothideomycetes_Taphrinomycotina_Chaetothyriales_Onygenales_18S <- c(Sordariomycetes_subset18S, Dothideomycetes_subset18S, Taphrinomycotina_subset18S, Chaetothyriales_subset18S, Onygenales_subset18S)
# Sordariomycetes_Dothideomycetes_Taphrinomycotina_Chaetothyriales_Onygenales_18S.msa <- msa(Sordariomycetes_Dothideomycetes_Taphrinomycotina_Chaetothyriales_Onygenales_18S, method = "ClustalOmega")
# Sordariomycetes_Dothideomycetes_Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay <- findAssay("Sordariomycetes_Dothideomycetes_Taphrinomycotina_Chaetothyriales_Onygenales", Sordariomycetes_Dothideomycetes_Taphrinomycotina_Chaetothyriales_Onygenales_18S.msa)
# #no assays


panfungal1 <- rbind(Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay[1,], Saccharomycetales_Basidiomycota_28S.assay[1,], Sordariomycetes_Dothideomycetes_18S.assay[1,], Mucorales_Eurotiales_18S.assay[1,])
write.csv(panfungal1, "panfungal_v1.csv", row.names = F)


# filtering. Hu cross reactivity
Hs18S <- Biostrings::readDNAStringSet("../BloodBorneFungalPathogens/assets/Hs18SrRNA.fas")

#Taphrinomycotina and Chaetothyriales and Onygenales
PrimerFor.cr <- unlist(sapply(unique(Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay$iupacSequenceFwd), function(p) vmatchPattern(p, Hs18S, max.mismatch = 3)))
PrimerFor.cr <- sapply(1:length(PrimerFor.cr), function(i) sapply(1:2, function(j) sum(PrimerFor.cr[[i]][[j]]@start)))[1,]

PrimerRev.cr <- unlist(sapply(unique(Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay$iupacSequenceRev), function(p) vmatchPattern(reverseComplement(DNAString(p)), Hs18S, max.mismatch = 3)))
PrimerRev.cr <- sapply(1:length(PrimerRev.cr), function(i) sapply(1:2, function(j) sum(PrimerRev.cr[[i]][[j]]@start)))[1,]

probe.cr <- unlist(sapply(unique(Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay$iupacSequencePr), function(p) vmatchPattern(p, Hs18S, max.mismatch = 3)))
probe.cr <- sapply(1:length(probe.cr), function(i) sapply(1:2, function(j) sum(probe.cr[[i]][[j]]@start)))[1,]

# selecting higher probe Tm
#Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay2 <- Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay2[Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay2$]
Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay2 <- as.data.frame(Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay) %>% dplyr::group_by(iupacSequenceFwd) %>% dplyr::filter(tmPr== max(unlist(tmPr))) %>% distinct

Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay3 <- Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay2[Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay2$tmPr > 60,]

# Saccharomycetales and Basidiomycota
Saccharomycetales_Basidiomycota_28S.assay 
PrimerFor.cr <- unlist(sapply(unique(Saccharomycetales_Basidiomycota_28S.assay$iupacSequenceFwd), function(p) vmatchPattern(p, Hs18S, max.mismatch = 3)))
PrimerFor.cr <- sapply(1:length(PrimerFor.cr), function(i) sapply(1:2, function(j) sum(PrimerFor.cr[[i]][[j]]@start)))[1,]

PrimerRev.cr <- unlist(sapply(unique(Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay$iupacSequenceRev), function(p) vmatchPattern(reverseComplement(DNAString(p)), Hs18S, max.mismatch = 3)))
PrimerRev.cr <- sapply(1:length(PrimerRev.cr), function(i) sapply(1:2, function(j) sum(PrimerRev.cr[[i]][[j]]@start)))[1,]

probe.cr <- unlist(sapply(unique(Saccharomycetales_Basidiomycota_28S.assay$iupacSequencePr), function(p) vmatchPattern(p, Hs18S, max.mismatch = 3)))
probe.cr <- sapply(1:length(probe.cr), function(i) sapply(1:2, function(j) sum(probe.cr[[i]][[j]]@start)))[1,]

Saccharomycetales_Basidiomycota_28S.assay$tmPr <- sapply(1:180, function(i) Saccharomycetales_Basidiomycota_28S.assay$tmPr[[i]][1])


tmp <- Saccharomycetales_Basidiomycota_28S.assay[1:2,]

for (i in 1:2){
  tmp$gcContentFwd[i] <- max(tmp$gcContentFwd[[i]])
  tmp$gcContentFwd[i] <- as.numeric(tmp$gcContentFwd[[i]])
}


Saccharomycetales_Basidiomycota_28S.assay2 <- as.data.frame(Saccharomycetales_Basidiomycota_28S.assay) %>% dplyr::group_by(iupacSequenceFwd) %>% dplyr::filter(tmPr== max(unlist(tmPr[1]))) %>% distinct
Saccharomycetales_Basidiomycota_28S.assay3 <- Saccharomycetales_Basidiomycota_28S.assay[Saccharomycetales_Basidiomycota_28S.assay2$tmPr > 56,]


# gblocks
TCO <- Taphrinomycotina_Chaetothyriales_Onygenales_18S.assay[1,]

Taphrinomycotina_Chaetothyriales_Onygenales_18S
TCOforMapping <- vmatchPattern(TCO$iupacSequenceFwd, Taphrinomycotina_Chaetothyriales_Onygenales_18S, max.mismatch = 1)
TCOfRevMapping <- vmatchPattern(reverseComplement(DNAString(TCO$iupacSequenceRev)), Taphrinomycotina_Chaetothyriales_Onygenales_18S, max.mismatch = 1)
TCOAmpliconLength <- unlist(sapply(1:57, function(i) TCOfRevMapping[[i]]@start - TCOforMapping[[i]]@start))
gblockTCO <- Taphrinomycotina_Chaetothyriales_Onygenales_18S[[1]][(TCOforMapping[[1]]@start-5):(TCOfRevMapping[[1]]@start + TCOfRevMapping[[1]]@width + 5)]

pairwiseAlignment(TCO$iupacSequenceFwd, gblockTCO, type = "local")
pairwiseAlignment(reverseComplement(DNAString(TCO$iupacSequenceRev)), gblockTCO, type = "local")
pairwiseAlignment(TCO$iupacSequencePr, gblockTCO, type = "local")

SB <- Saccharomycetales_Basidiomycota_28S.assay[1,]
SBforMapping <- vmatchPattern(SB$iupacSequenceFwd, Saccharomycetales_Basidiomycota_28S, max.mismatch = 1)
SBrevMapping <- vmatchPattern(reverseComplement(DNAString(SB$iupacSequenceRev)), Saccharomycetales_Basidiomycota_28S, max.mismatch = 1)
SBAmpliconLength <- unlist(sapply(1:277, function(i) SBrevMapping[[i]]@start - SBforMapping[[i]]@start))

gblockSB <- Saccharomycetales_Basidiomycota_28S[[1]][(SBforMapping[[1]]@start-5):(SBrevMapping[[1]]@start + SBrevMapping[[1]]@width + 5)]

pairwiseAlignment(SB$iupacSequenceFwd, gblockSB, type = "local")
pairwiseAlignment(reverseComplement(DNAString(SB$iupacSequenceRev)), gblockSB, type = "local")
pairwiseAlignment(SB$iupacSequencePr, gblockSB, type = "local")


SD <- Sordariomycetes_Dothideomycetes_18S.assay[1,]
SDforMapping <- vmatchPattern(SD$iupacSequenceFwd, Sordariomycetes_Dothideomycetes_18S, max.mismatch = 0)
SDrevMapping <- vmatchPattern(reverseComplement(DNAString(SD$iupacSequenceRev)), Sordariomycetes_Dothideomycetes_18S, max.mismatch = 0)
SBAmpliconLength <- unlist(sapply(1:176, function(i) SDrevMapping[[i]]@start - SDforMapping[[i]]@start))
gblockSD <- Sordariomycetes_Dothideomycetes_18S[[1]][(SDforMapping[[1]]@start-5):(SDrevMapping[[1]]@start + SDrevMapping[[1]]@width + 5)]

pairwiseAlignment(SD$iupacSequenceFwd, gblockSD, type = "local")
pairwiseAlignment(reverseComplement(DNAString(SD$iupacSequenceRev)), gblockSD, type = "local")
pairwiseAlignment(SD$iupacSequencePr, gblockSD, type = "local")


ME <- Mucorales_Eurotiales_18S.assay[1,]
MEforMapping <- vmatchPattern(ME$iupacSequenceFwd, Mucorales_Eurotiales_18S, max.mismatch = 0)
MErevMapping <- vmatchPattern(reverseComplement(DNAString(ME$iupacSequenceRev)), Mucorales_Eurotiales_18S, max.mismatch = 0)
MEAmpliconLength <- unlist(sapply(1:268, function(i) MErevMapping[[i]]@start - MEforMapping[[i]]@start))
gblockME <- Mucorales_Eurotiales_18S[[1]][(MEforMapping[[1]]@start-5):(MErevMapping[[1]]@start + MErevMapping[[1]]@width + 5)]

pairwiseAlignment(ME$iupacSequenceFwd, gblockME, type = "local")
pairwiseAlignment(reverseComplement(DNAString(ME$iupacSequenceRev)), gblockME, type = "local")
pairwiseAlignment(ME$iupacSequencePr, gblockME, type = "local")

PanFungalgBlocks <- DNAStringSet(list(gblockTCO, gblockSB, gblockSD, gblockME))
names(PanFungalgBlocks) <- c("Taphrinomycotina_Chaetothyriales_Onygenales", "Saccharomycetales_Basidiomycota", "Sordariomycetes_Dothideomycetes", "Mucorales_Eurotiales")
writeXStringSet(PanFungalgBlocks, "PanFungalgBlocks.fas")


### April 26, 2023: double checking specificity of pan-Fungal reagents###
Hs18S <- Biostrings::readDNAStringSet("../BloodBorneFungalPathogens/assets/Hs18SrRNA.fas")
# filtering. Hu cross reactivity
Hs28S <- Biostrings::readDNAStringSet("../BloodBorneFungalPathogens/assets/Hs28SrRNA.fas")

pan_fugnal_v1 <- read.csv("panfungal_v1.csv")

# Taphrinomycotina_Chaetothyriales_Onygenales (18S)
TCO_forPrimer.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Taphrinomycotina_Chaetothyriales_Onygenales",]$iupacSequenceFwd, Hs18S, max.mismatch = 3)
TCO_forPrimer.ali <- pairwiseAlignment(pan_fugnal_v1[pan_fugnal_v1$fungi == "Taphrinomycotina_Chaetothyriales_Onygenales",]$iupacSequenceFwd, Hs18S[[1]], type = "local")

TCO_probe.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Taphrinomycotina_Chaetothyriales_Onygenales",]$iupacSequencePr, Hs18S, max.mismatch = 3)

TCO_revPrimer.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Taphrinomycotina_Chaetothyriales_Onygenales",]$iupacSequenceRev, reverseComplement(Hs18S), max.mismatch = 3)

# Saccharomycetales_Basidiomycota (28S)
SB_forPrimer.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Saccharomycetales_Basidiomycota",]$iupacSequenceFwd, Hs28S, max.mismatch = 3)

SB_probe.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Saccharomycetales_Basidiomycota",]$iupacSequencePr, Hs28S, max.mismatch = 3)

SB_revPrimer.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Saccharomycetales_Basidiomycota",]$iupacSequenceRev, reverseComplement(Hs28S), max.mismatch = 3)

# Sordariomycetes_Dothideomycetes (18S)
SD_forPrimer.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Sordariomycetes_Dothideomycetes",]$iupacSequenceFwd, Hs18S, max.mismatch = 1)

SD_probe.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Sordariomycetes_Dothideomycetes",]$iupacSequencePr, Hs18S, max.mismatch = 4)

SD_revPrimer.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Sordariomycetes_Dothideomycetes",]$iupacSequenceRev, reverseComplement(Hs18S), max.mismatch = 4)

# Mucorales_Eurotiales (18S)
ME_forPrimer.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Mucorales_Eurotiales",]$iupacSequenceFwd, Hs18S, max.mismatch = 3)

ME_probe.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Mucorales_Eurotiales",]$iupacSequenceRev, Hs18S, max.mismatch = 3)

ME_probe2.cr <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Mucorales_Eurotiales",]$iupacSequenceRev, Hs28S, max.mismatch = 3)


## April 26, 2023: double checking inclusivity of pan-Fungal reagents###
TCO_forPrimer.map <- vmatchPattern(pan_fugnal_v1[pan_fugnal_v1$fungi == "Taphrinomycotina_Chaetothyriales_Onygenales",]$iupacSequenceFwd, Taphrinomycotina_Chaetothyriales_Onygenales_18S, max.mismatch = 1)



