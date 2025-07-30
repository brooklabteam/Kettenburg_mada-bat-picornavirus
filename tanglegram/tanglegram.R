rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)
library(phylotools)
library(phytools)
library(phylobase)
library(cowplot)
library(ggtreeExtra)
library(viridis)
library(ggplotify)
library(patchwork)

###packages loaded
##########################################################################################################
##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_mada-bat-picornavirus/"
setwd(paste0(homewd,"/tanglegram"))

tr1 <- read.tree("bat_sapelo_hosts_align_fasttree.newick")
tr2 <- read.tree("batsapelo_fullname_align_fasttree.newick")

assoc<-cbind(c("E_helvum","E_dupreanum", "E_dupreanum","R_madagascariensis","E_spelaea",
               "E_spelaea","E_spelaea","R_leschenaultii","R_aegyptiacus","E_helvum","E_helvum"),
             c("NC_033820.1","OQ818320","OQ818321","OQ818329","OR951325.1","OR951326.1",
               "OR951327.1","OR951332.1","PP711911.1","PP711921.1","PP711943.1"))

obj<-cophylo(tr1,tr2,assoc,rotate=FALSE)
plot(obj)


