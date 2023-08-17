rm(list=ls())

library(plyr)
library(dplyr)
library(seqinr)


#set wd
homewd = "/Users/gwenddolenkettenburg/developer/mada-bat-picornavirus"
setwd(paste0(homewd, "/BEAST/"))

#load the dataset and query
dat <- read.csv(file = "picornaviridae_beast_metadata.csv", header = T, stringsAsFactors = F)
head(dat)

#now, get the fasta file 
fasta.dat <- read.fasta("compiled_seq.fasta", forceDNAtolower = F, as.string = T)

names(fasta.dat)

fasta.meta <-  cbind.data.frame(tip_label = names(fasta.dat))

#add beast name to the main dataset
dat$Collection_Date <- as.Date(dat$Collection_Date)
#dat$collection_date <- as.character(dat$collection_date)
#dat$collection_date[is.na(dat$collection_date)] <- paste0(dat$collection_year[is.na(dat$collection_date)], "07-31")
#dat$collection_date <- as.Date(dat$collection_date)

dat$beast_name <- paste0(dat$Accession, "_", dat$Collection_Date)

dat.merge <- dplyr::select(dat, tip_label, beast_name)


setdiff(dat.merge$tip_label, fasta.meta$tip_label)
setdiff(fasta.meta$tip_label,dat.merge$tip_label)

#fasta.meta$tip_label[fasta.meta$tip_label=="DQ648794_1_Bat_coronavirus_(BtCoV_133_2005)"] <- "DQ648794_1_Bat_coronavirus__BtCoV_133_2005"

#and add to data
fasta.meta <- merge(fasta.meta, dat.merge, all.x = T, by="tip_label", sort = F)
head(fasta.meta)
fasta.meta$beast_name
#and write over:
write.fasta(fasta.dat, names=fasta.meta$beast_name, "picornaviridae_beast.fasta")

