
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
library(devtools)
library(tangler)

###packages loaded
##########################################################################################################
##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_mada-bat-picornavirus/"
setwd(paste0(homewd,"/tanglegram"))

#bat picornavirus
#load the tree and root it
t1 <-  read.tree("batpicorna_host_nexfile") 
t2 <-  read.tree("batpicorna_virus_nexfile") 

#load tree data prepared from elsewhere
dat <- read.csv(("bat_picornaviruses.csv"), header = T, stringsAsFactors = F)
head(dat)

dat$event<-factor(dat$event)

# Annotate trees
tree1 <- ggtree(t1)   %<+% dat +
  geom_tiplab() +
  geom_tippoint(aes(color=event))
tree1


tree2 <- ggtree(t2) %<+% dat + geom_tiplab()
tree2

simple.tanglegram(tree1, tree2, event, A, tiplab = T)

simple.tanglegram(tree1, tree2, event, A, l_color = 'green3',  t2_pad = 0.3,
                  tiplab = T, lab_pad = 0.1, x_hjust = 1, t2_tiplab_size = 3)


