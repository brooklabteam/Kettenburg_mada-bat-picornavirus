rm(list=ls())

library(plyr)
library(dplyr)


#set wd
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus"
setwd(paste0(homewd, "/fasta_files/NCBI/reference_sequences/picornaviridae"))

#load the dataset and query - change dataset and working directory as you go
dat <- read.csv(file = "teschovirus_fullgenome_nt_ncbi.csv", header = T, stringsAsFactors = F)

head(dat)

#and get the text to download from NCBI, check duplicated accession numbers
dat <- dat[!duplicated(dat$Accession),] #no duplicates
accession_num <- paste(c(dat$Accession), collapse = ",")

#now put this into your webbrowser to download
text.for.NCBI <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=",accession_num)
text.for.NCBI

###YOU ONLY NEED TO DO THE BELOW STEP IF YOU DO NOT ALIGN IN GENEIOUS AND INSTEAD USE THE MAFFT WEBSITE

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
