rm(list=ls())

library(plyr)
library(dplyr)


#set wd
homewd = "/Users/gwenddolenkettenburg/Desktop/mada_bat_picornavirus"
setwd(paste0(homewd, "/fasta_files/NCBI/"))

#load the dataset and query
dat <- read.csv(file = "picornaviridae_refseq_metadata.csv", header = T, stringsAsFactors = F)

head(dat)
#look at the unique hosts
sort(unique(dat$Host))
sort(unique(dat$Species[dat$Host==""])) #these ones have no hosts

#select all the bats
bat.dat = subset(dat, Host=="Chiroptera" | Host =="Hipposideros armiger" | Host =="Miniopterus magnater"|
                   Host=="Miniopterus pusillus" | Host == "Miniopterus schreibersii" | Host =="Myotis ricketti" | Host =="Pipistrellus pipistrellus"|
                   Host== "Rhinolophus sinicus" | Host== "Eidolon helvum" | Host== "Miniopterus" | Host=="Miniopterus fuliginosus")

sort(unique(bat.dat$Species))

nrow(bat.dat) 

#then select the ref seq for the hosts that are not bats
ref.dat <- subset(dat, Host=="") #check those with no host, have to manually fill in 

#now take just those with no bats...
ref.sub = subset(dat, Host!="Chiroptera" & Host!="Hipposideros armiger" & Host !="Miniopterus magnater" & Host !="Miniopterus pusillus" &
                   Host!="Miniopterus schreibersii" & Host!="Myotis ricketti" & Host !="Pipistrellus pipistrellus" & Host !="Rhinolophus sinicus" &
                   Host!="Eidolon helvum" & Host!="Miniopterus" & Host!="Miniopterus fuliginosus")

#now put these together to draw from GenBank
bat.picorna <- rbind(bat.dat, ref.sub)

#Check duplicated records
bat.picorna <- bat.picorna[!duplicated(bat.picorna),] #no duplicates

#and get the text to download from NCBI, check duplicated accession numbers
bat.picorna <- bat.picorna[!duplicated(bat.picorna$Accession),] #no duplicates


#and remove those that are repeats of the same record 
#all.bat.picorna = subset(all.CoV, Accession!="KU762338" &
                   # Accession!= "KF636752" &
                   # Accession!="GU190215" &
                   # Accession!="EF065505" &
                   # Accession!="EF065509"&
                   # Accession!="EF065513" &
                   # Accession!="KX574227")

accession_num <- paste(c(bat.picorna$Accession), collapse = ",")

#now put this into your webbrowser to download
text.for.NCBI <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=",accession_num)

#once downloaded, send to MAFFT for alignment
rm(list=ls())
#then, after alignment is ready, prepare the names for RAxML (no space, semicolon, colon, parentheses, dash, slash, comma, quote allowed in name (should just all be underscore)
library(seqinr)
#library(msa)
alignment1 <- read.alignment(file = "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/fasta_files/NCBI/trim_align_picornaviridae.fasta", format="fasta", forceToLower = F)


tmp <- as.list(alignment1$nam)

change.spacing <- function(df){
  df_new <- sapply(strsplit(df,"-"), function(x) x[[1]])
  return(df_new)
}

names_new = c(unlist(lapply(tmp, change.spacing)))

#new_names <- sub("__", "_", new_names) 
class(alignment1$seq)
write.fasta(sequences = as.list(alignment1$seq), names = names_new, file.out =  "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/fasta_files/NCBI/cleaned_picornaviridae.fasta", as.string = T, open="w")

#now send to modeltest and eventually RAxML
