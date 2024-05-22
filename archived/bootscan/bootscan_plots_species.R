rm(list=ls())

library(ggplot2)
library(gggenes)
library(cowplot)
library(gridExtra)
library(grid)
library(plyr)
library(dplyr)
library(reshape2)
library(cowplot)
library(patchwork)
library(ggplotify)


#This is to make a figure of the representative bootscan plots with their matching genome plots

##First make the genome plots
homewd="/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/genome_annotation_and_characterization/"
setwd("~/Desktop/developer/mada-bat-picornavirus/genome_annotation_and_characterization/genus_gene_maps")

#Load the gene data
ictv <- read.csv("ictv_blast_genes_species.csv", header = T, stringsAsFactors = F)
ictv$gene<-factor(ictv$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                        "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                        "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Minor structural protein", "Non-structural polyprotein",
                                        "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                        "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the peptide files
ictv_pep <- read.csv("ictv_blast_peptides_species.csv", header = T, stringsAsFactors = F)
ictv_pep$gene<-factor(ictv_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                                "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Minor structural protein", "Non-structural polyprotein",
                                                "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the feature file in case its needed
ictv_feat <- read.csv("ictv_blast_features_species.csv", header = T, stringsAsFactors = F)
ictv_feat$gene<-factor(ictv_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                  "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                                  "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Minor structural protein", "Non-structural polyprotein",
                                                  "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                  "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))

#Pick colors for genes
colz=c("5'UTR"="gold", "L"="royalblue","VP4"="paleturquoise3", "VP2"="skyblue1", "VP0"="royalblue4", "VP3"="steelblue1",
       "VP1"="cadetblue1", "VP1/2A"="cadetblue1", "2A"="orange1", "2B"="sienna2", "2C"="darkorange1", "3A"="palevioletred1", "3B"="plum",
       "3C"="rosybrown1", "3D"="pink2", "Helicase"="darkseagreen1","NS4"="darkolivegreen1","Vpg"="seagreen1","Pro-Pol"="palegreen2",
       "Polyprotein"="azure3","Putative polyprotein"="mediumorchid1", "Non-structural polyprotein"="mediumorchid4", 
       "Minor structural protein" ="slateblue3", "Structural polyprotein"="lightgoldenrod1",
       "Hypothetical protein"="darkslategrey", "Similar to structural polyprotein"="mediumpurple1", "Similar to putative polyprotein"="mediumpurple3", 
       "Similar to polyprotein"="mediumpurple4","3'UTR"="yellow")


#Plot ICTV and BLAST together plots

#Subset ICTV data by species
ictv_eidolon<-subset(ictv,molecule=="Eidolon")
ictv_eidolon_full<-subset(ictv_eidolon,type=="full")
ictv_eidolon_p1<-subset(ictv_eidolon,type=="p1")
ictv_pteropus<-subset(ictv,molecule=="Pteropus")
ictv_rousettus<-subset(ictv,molecule=="Rousettus")
ictv_rousettus_full<-subset(ictv_rousettus,type=="full")
ictv_rousettus_all<-subset(ictv_rousettus,type=="all")

ictv_eidolon_pep<-subset(ictv_pep,molecule=="Eidolon")
ictv_eidolon_full_pep<-subset(ictv_eidolon_pep,type=="full")
ictv_eidolon_p1_pep<-subset(ictv_eidolon_pep,type=="p1")
ictv_pteropus_pep<-subset(ictv_pep,molecule=="Pteropus")
ictv_rousettus_pep<-subset(ictv_pep,molecule=="Rousettus")
ictv_rousettus_full_pep<-subset(ictv_rousettus_pep,type=="full")
ictv_rousettus_all_pep<-subset(ictv_rousettus_pep,type=="all")

ictv_eidolon_feat<-subset(ictv_feat,molecule=="Eidolon")
ictv_eidolon_full_feat<-subset(ictv_eidolon_feat,type=="full")
ictv_eidolon_p1_feat<-subset(ictv_eidolon_feat,type=="p1")
ictv_pteropus_feat<-subset(ictv_feat,molecule=="Pteropus")
ictv_rousettus_feat<-subset(ictv_feat,molecule=="Rousettus")
ictv_rousettus_full_feat<-subset(ictv_rousettus_feat,type=="full")
ictv_rousettus_all_feat<-subset(ictv_rousettus_feat,type=="all")


#plot ictv and blast plots
ictv_eidolon_full<-ggplot(ictv_eidolon_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_eidolon_full_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                  xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_eidolon_full_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,10600),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_eidolon_full


ictv_eidolon_p1<-ggplot(ictv_eidolon_p1, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_eidolon_p1_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                               xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_eidolon_p1_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1190,2740),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_eidolon_p1


ictv_pteropus<-ggplot(ictv_pteropus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_pteropus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                               xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_pteropus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1,5180),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_pteropus


ictv_rousettus_all<-ggplot(ictv_rousettus_all, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_rousettus_all_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                           xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_rousettus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(5360,6500),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_rousettus_all


ictv_rousettus_full<-ggplot(ictv_rousettus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_rousettus_full_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_rousettus_full_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1,8900),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_rousettus_full



##Now get the plots from bootscan

setwd("~/Desktop/developer/mada-bat-picornavirus/bootscan/species_specific_output")

colzpalette<-c("#8ECAE6","#E48B97","#219EBC","#B52B09","#023047","#A60067","#FFB703","#987B6F","#FB8500","#FB8500", "orchid")


#Eidolon full
eidolon_full_bootscan <- read.csv(file = "eidolon_full_bootscan.csv", header = T, stringsAsFactors = F)
head(eidolon_full_bootscan)

#move to long
long.sim_nt <- melt(eidolon_full_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818317","OQ818318","OQ818321","OQ818337", "OQ818319", "OQ818322"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818317"] <- "Eidolon dupreanum kunsagivirus OQ818317"
long.sim_nt$accession[long.sim_nt$accession == "OQ818318"] <- "Eidolon dupreanum teschovirus OQ818318"
long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_nt$accession[long.sim_nt$accession == "OQ818337"] <- "Eidolon dupreanum hepatovirus OQ818337"
long.sim_nt$accession[long.sim_nt$accession == "OQ818319"] <- "Eidolon dupreanum sapovirus OQ818319"
long.sim_nt$accession[long.sim_nt$accession == "OQ818322"] <- "Eidolon dupreanum kobuvirus OQ818322"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum kunsagivirus OQ818317", 
                                                                  "Eidolon dupreanum teschovirus OQ818318",
                                                                  "Eidolon dupreanum sapelovirus OQ818321",
                                                                  "Eidolon dupreanum hepatovirus OQ818337",
                                                                  "Eidolon dupreanum sapovirus OQ818319",
                                                                  "Eidolon dupreanum kobuvirus OQ818322"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapelovirus "), "OQ818320"))

eidolon_full <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 4))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
  #                    labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

eidolon_full

#put gene map with PySimPlot
edup_full<-eidolon_full/ictv_eidolon_full+plot_layout(nrow=2,  heights = c(2, 0.30))
edup_full

edup_full<-as.ggplot(edup_full)
edup_full



#Eidolon p1
eidolon_p1_bootscan <- read.csv(file = "eidolon_p1_bootscan.csv", header = T, stringsAsFactors = F)
head(eidolon_p1_bootscan)

#move to long
long.sim_nt <- melt(eidolon_p1_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818317","OQ818318","OQ818321","OQ818337",
                                                                                  "OQ818319", "OQ818322", "OQ818340","OQ818343",
                                                                                  "OQ818339","OQ818341","OQ818338"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818317"] <- "Eidolon dupreanum kunsagivirus OQ818317"
long.sim_nt$accession[long.sim_nt$accession == "OQ818318"] <- "Eidolon dupreanum teschovirus OQ818318"
long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_nt$accession[long.sim_nt$accession == "OQ818337"] <- "Eidolon dupreanum hepatovirus OQ818337"
long.sim_nt$accession[long.sim_nt$accession == "OQ818319"] <- "Eidolon dupreanum sapovirus OQ818319"
long.sim_nt$accession[long.sim_nt$accession == "OQ818322"] <- "Eidolon dupreanum kobuvirus OQ818322"
long.sim_nt$accession[long.sim_nt$accession == "OQ818340"] <- "Eidolon dupreanum sapovirus OQ818340"
long.sim_nt$accession[long.sim_nt$accession == "OQ818343"] <- "Eidolon dupreanum sapelovirus OQ818343"
long.sim_nt$accession[long.sim_nt$accession == "OQ818339"] <- "Eidolon dupreanum picorna-like virus OQ818339"
long.sim_nt$accession[long.sim_nt$accession == "OQ818341"] <- "Eidolon dupreanum felisavirus OQ818341"
long.sim_nt$accession[long.sim_nt$accession == "OQ818338"] <- "Eidolon dupreanum cripavirus OQ818338"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum kunsagivirus OQ818317", 
                                                                  "Eidolon dupreanum teschovirus OQ818318",
                                                                  "Eidolon dupreanum sapelovirus OQ818321",
                                                                  "Eidolon dupreanum hepatovirus OQ818337",
                                                                  "Eidolon dupreanum sapovirus OQ818319",
                                                                  "Eidolon dupreanum kobuvirus OQ818322", 
                                                                  "Eidolon dupreanum sapovirus OQ818340",
                                                                  "Eidolon dupreanum sapelovirus OQ818343",
                                                                  "Eidolon dupreanum picorna-like virus OQ818339",
                                                                  "Eidolon dupreanum felisavirus OQ818341",
                                                                  "Eidolon dupreanum cripavirus OQ818338"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapelovirus "), "OQ818320"))

eidolon_p1 <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 4))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
  #                    labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

eidolon_p1

#put gene map with PySimPlot
edup_p1<-eidolon_p1/ictv_eidolon_p1+plot_layout(nrow=2,  heights = c(2, 0.30))
edup_p1

edup_p1<-as.ggplot(edup_p1)
edup_p1


#Pteropus
pteropus_bootscan <- read.csv(file = "pteropus_all_bootscan.csv", header = T, stringsAsFactors = F)
head(pteropus_bootscan)

#move to long
long.sim_nt <- melt(pteropus_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818316","OQ818332","OQ818336","OQ818334",
                                                                                  "OQ818330", "OQ818335", "OQ818333"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818316"] <- "Pteropus rufus mischivirus OQ818316"
long.sim_nt$accession[long.sim_nt$accession == "OQ818332"] <- "Pteropus rufus bee virus OQ818332"
long.sim_nt$accession[long.sim_nt$accession == "OQ818336"] <- "Pteropus rufus picorna-like virus OQ818336"
long.sim_nt$accession[long.sim_nt$accession == "OQ818334"] <- "Pteropus rufus picorna-like virus OQ818334"
long.sim_nt$accession[long.sim_nt$accession == "OQ818330"] <- "Pteropus rufus cheravirus OQ818330"
long.sim_nt$accession[long.sim_nt$accession == "OQ818335"] <- "Pteropus rufus felisavirus OQ818335"
long.sim_nt$accession[long.sim_nt$accession == "OQ818333"] <- "Pteropus rufus picorna-like virus OQ818333"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Pteropus rufus mischivirus OQ818316", 
                                                                  "Pteropus rufus bee virus OQ818332",
                                                                  "Pteropus rufus picorna-like virus OQ818336",
                                                                  "Pteropus rufus picorna-like virus OQ818334",
                                                                  "Pteropus rufus cheravirus OQ818330",
                                                                  "Pteropus rufus felisavirus OQ818335", 
                                                                  "Pteropus rufus picorna-like virus OQ818333"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("Reference: ",italic("Pteropus rufus tetnovirus "), "OQ818331"))

pteropus <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 4))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
  #                    labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

pteropus

#put gene map with PySimPlot
pruf<-pteropus/ictv_pteropus+plot_layout(nrow=2,  heights = c(2, 0.30))
pruf

pruf<-as.ggplot(pruf)
pruf



#rousettus full
rousettus_full_bootscan <- read.csv(file = "rousettus_full_bootscan.csv", header = T, stringsAsFactors = F)
head(rousettus_full_bootscan)

#move to long
long.sim_nt <- melt(rousettus_full_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818323","OQ818324","OQ818325","OQ818329"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818323"] <- "Rousettus madagascariensis teschovirus OQ818323"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324"] <- "Rousettus madagascariensis teschovirus OQ818324"
long.sim_nt$accession[long.sim_nt$accession == "OQ818325"] <- "Rousettus madagascariensis picornavirus OQ818325"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rousettus madagascariensis teschovirus OQ818323", 
                                                                  "Rousettus madagascariensis teschovirus OQ818324",
                                                                  "Rousettus madagascariensis picornavirus OQ818325",
                                                                  "Rousettus madagascariensis sapelovirus OQ818329"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("Reference: ",italic("Rousettus madagascariensis picornavirus "), "OQ818328"))

rousettus_full <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 4))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
  #                    labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

rousettus_full

#put gene map with PySimPlot
rou_full<-rousettus_full/ictv_rousettus_full+plot_layout(nrow=2,  heights = c(2, 0.30))
rou_full

rou_full<-as.ggplot(rou_full)
rou_full



#rousettus all
rousettus_all_bootscan <- read.csv(file = "rousettus_all_bootscan.csv", header = T, stringsAsFactors = F)
head(rousettus_all_bootscan)

#move to long
long.sim_nt <- melt(rousettus_all_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818323","OQ818324","OQ818325","OQ818328",
                                                                                  "OQ818346", "OQ818329", "OQ818345","OQ818347"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818323"] <- "Rousettus madagascariensis teschovirus OQ818323"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324"] <- "Rousettus madagascariensis teschovirus OQ818324"
long.sim_nt$accession[long.sim_nt$accession == "OQ818325"] <- "Rousettus madagascariensis picornavirus OQ818325"
long.sim_nt$accession[long.sim_nt$accession == "OQ818328"] <- "Rousettus madagascariensis picornavirus OQ818328"
long.sim_nt$accession[long.sim_nt$accession == "OQ818346"] <- "Rousettus madagascariensis picornavirus OQ818346"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_nt$accession[long.sim_nt$accession == "OQ818345"] <- "Rousettus madagascariensis sapovirus OQ818345"
long.sim_nt$accession[long.sim_nt$accession == "OQ818347"] <- "Rousettus madagascariensis sapovirus OQ818347"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rousettus madagascariensis teschovirus OQ818323", 
                                                                  "Rousettus madagascariensis teschovirus OQ818324",
                                                                  "Rousettus madagascariensis picornavirus OQ818325",
                                                                  "Rousettus madagascariensis picornavirus OQ818328",
                                                                  "Rousettus madagascariensis picornavirus OQ818346",
                                                                  "Rousettus madagascariensis sapelovirus OQ818329", 
                                                                  "Rousettus madagascariensis sapovirus OQ818345",
                                                                  "Rousettus madagascariensis sapovirus OQ818347"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("Reference: ",italic("Rousettus madagascariensis sapovirus "), "OQ818348"))

rousettus_all <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 4))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
  #                    labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

rousettus_all

#put gene map with PySimPlot
rou_all<-rousettus_all/ictv_rousettus_all+plot_layout(nrow=2,  heights = c(2, 0.30))
rou_all

rou_all<-as.ggplot(rou_all)
rou_all


##Now put the whole figure together

#all bootscan together
fig<-plot_grid(edup_full,
                           edup_p1,
                           pruf,
                           rou_all,
                    ncol=2,
                    labels="AUTO",  label_size = 23, align = "hv", axis="b")
fig

#export landscape 15x20 inch PDF

