rm(list=ls())


library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(lubridate)
#library(rBt)
library(treeio)


#make Bayesian timetree from caliciviridae strict molecular clock model

#first, read in the tree

homewd= "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/"

setwd(paste0(homewd, "/BEAST/refseq_caliciviridae/strict"))

tree <-  read.beast(file = paste0(homewd, "/BEAST/refseq_caliciviridae/strict/caliciviridae_strict_mean_tree"))

treedat <- cbind.data.frame(tip_name = tree@phylo$tip.label)
treedat$beast_name <-treedat$tip_name
#tree <- read.annot.beast(file = paste0(homewd, "/Fig4/beast-out/AllNobeco/NobecoStrict/AvgNobecoStrictNexus.trees"))

#tree$node.label <- round(tree$posterior,2)
#treedat <- cbind.data.frame(tip_name = tree$tip.label)
treedat$accession_num <- sapply(strsplit(treedat$tip_name, "_"), function(x) x[[1]])
treedat$accession_num[treedat$accession_num=="NC"] <- c("NC_000940","NC_001481","NC_001543","NC_001959","NC_002551","NC_002615","NC_004064",
                                                        "NC_004541","NC_004542","NC_006269","NC_006554","NC_006875","NC_007916","NC_008311",
                                                        "NC_008580","NC_010624","NC_011050","NC_011704","NC_012699","NC_017936","NC_019712",
                                                        "NC_024031","NC_024078","NC_025676","NC_027026","NC_027122","NC_029645","NC_029646",
                                                        "NC_029647","NC_030793","NC_031324","NC_033081","NC_033776","NC_034444","NC_035675",
                                                        "NC_039475","NC_039476","NC_039477","NC_039897","NC_040674","NC_040876","NC_043512",
                                                        "NC_043516","NC_044045","NC_044046","NC_044047","NC_044853","NC_044854","NC_044855",
                                                        "NC_044856","NC_044932","NC_045762")
treedat$accession_num[treedat$accession_num=="F"] <- c("OQ818319")

#names(treedat)[names(treedat)=="tip_name"] <- "beast_name"

#and load data of corresponding tree

dat <- read.csv(file = "caliciviridae_beast_metadata.csv", header = T, stringsAsFactors = F)
dat$Collection_Date <- as.Date(dat$Collection_Date)


colz = c("Sapovirus" = "royalblue3",    "Vesivirus"  = "turquoise1",   "Lagovirus"  = "goldenrod1",   "Norovirus"   = "dodgerblue1" ,   "Calicivirus" = "firebrick1" ,
         "Salovirus"  = "lightpink1" ,    "Bavovirus"  = "hotpink1" ,  "Minovirus" = "lightskyblue" ,   "Coronavirus"  = "black", "Recovirus"  = "darkorange1", 
         "Nacovirus"="thistle3", "Nebovirus"="darkorchid4")

#pick order for the labels
dat$Genus <- factor(dat$Genus, levels = c("Sapovirus" ,  "Vesivirus",   "Lagovirus",   "Norovirus",   "Calicivirus",
                                          "Salovirus",    "Bavovirus",  "Minovirus",  "Recovirus", "Nacovirus", "Nebovirus", 
                                          "Coronavirus"))   

dat$novel <- as.factor(dat$novel)


mrsd.dat <- max(dat$Collection_Date)
p1 <- ggtree(tree, mrsd=mrsd.dat)  + theme_tree2()  +geom_nodelab()
p1

tree.dat <- p1$data
node.sub <- dplyr::select(tree.dat, node, x)
names(node.sub) <-  c("node", "nodetime")

#and 
head(dat)

dat.plot <- merge(treedat, dat, by="accession_num", all.x = T, sort=F)

head(dat.plot)
dat.plot$new_label = NA
dat.plot$new_label[!is.na(dat.plot$Species)] <- paste(dat.plot$accession_num[!is.na(dat.plot$Species)], " | ", 
                                                     dat.plot$Species[!is.na(dat.plot$Species)], " | ", 
                                                     dat.plot$source[!is.na(dat.plot$Species)], " | ",
                                                     dat.plot$country[!is.na(dat.plot$Species)], " | ",
                                                     dat.plot$collection_year[!is.na(dat.plot$Species)])

dat.plot$new_label[is.na(dat.plot$Species)] <- paste(dat.plot$accession_num[is.na(dat.plot$Species)], " | ", 
                                                    dat.plot$source[is.na(dat.plot$Species)], " | ",
                                                    dat.plot$country[is.na(dat.plot$Species)], " | ",
                                                    dat.plot$collection_year[is.na(dat.plot$Species)])


tree@phylo$tip.label <- dat.plot$new_label

dat.sub <- dplyr::select(dat.plot, new_label, Genus, Family, Collection_Date,Host, Country, source, Collection_Year, novel, Isolate, Species)
head(dat.sub)

# dat.sub$novel = "no"
# dat.sub$novel[dat.sub$country=="Madagascar"] <- "yes"
# 
# colz2 = c('yes' =  "yellow", 'no' = "white")

dat.sub$Host[dat.sub$Host==0] <- "Non-bat host"
dat.sub$Host[dat.sub$Host==1] <- "Bat host"
dat.sub$Host <- as.factor(dat.sub$Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")

#plot tree
p1 <-ggtree(tree, mrsd=mrsd.dat, size=.8) %<+% dat.sub+
  geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +  
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  #theme_tree2() +
  coord_cartesian(clip = "off") + 
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=3, stroke=.1)+
  scale_fill_continuous(low="yellow", high="red") +
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1)) +
  theme(legend.position = "left", 
        legend.direction = "vertical",
        legend.text = element_text(size=8), 
        legend.key.size = unit(0.3, "cm"))+
  scale_x_continuous(breaks=c(1400, 1600, 1800, 2000),
                     labels=c(623, 423, 223, 23)) +
  xlab("Years to MRCA")

p1


