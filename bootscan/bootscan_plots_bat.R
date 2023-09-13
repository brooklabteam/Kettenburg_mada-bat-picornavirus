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
#Load the gene data
ictv <- read.csv("ictv_blast_genes.csv", header = T, stringsAsFactors = F)
ictv$gene<-factor(ictv$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                        "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                        "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                        "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                        "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the peptide files
ictv_pep <- read.csv("ictv_blast_peptides.csv", header = T, stringsAsFactors = F)
ictv_pep$gene<-factor(ictv_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                                "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                                "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the feature file in case its needed
ictv_feat <- read.csv("ictv_blast_features.csv", header = T, stringsAsFactors = F)
ictv_feat$gene<-factor(ictv_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                  "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                                  "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                                  "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                  "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))

#Pick colors for genes
colz=c("5'UTR"="gold", "L"="royalblue","VP4"="paleturquoise3", "VP2"="skyblue1", "VP0"="royalblue4", "VP3"="steelblue1",
       "VP1"="cadetblue1", "VP1/2A"="cadetblue1", "2A"="orange1", "2B"="sienna2", "2C"="darkorange1", "3A"="palevioletred1", "3B"="plum",
       "3C"="rosybrown1", "3D"="pink2", "Helicase"="darkseagreen1","NS4"="darkolivegreen1","Vpg"="seagreen1","Pro-Pol"="palegreen2",
       "Polyprotein"="azure3","Putative polyprotein"="mediumorchid1", "Non-structural polyprotein"="mediumorchid4", 
       "Putative minor structural protein" ="slateblue3", "Structural polyprotein"="lightgoldenrod1",
       "Hypothetical protein"="darkslategrey", "Similar to structural polyprotein"="mediumpurple1", "Similar to putative polyprotein"="mediumpurple3", 
       "Similar to polyprotein"="mediumpurple4","3'UTR"="yellow")


#Plot ICTV and BLAST together plots

#Subset ICTV data by virus
ictv_batpicornavirus<-subset(ictv,molecule=="Bat picornavirus")
ictv_hepatovirus<-subset(ictv,molecule=="Hepatovirus")
ictv_kobuvirus<-subset(ictv,molecule=="Kobuvirus")
ictv_kunsagivirus<-subset(ictv,molecule=="Kunsagivirus")
ictv_mischivirus<-subset(ictv,molecule=="Mischivirus")
ictv_sapovirus<-subset(ictv,molecule=="Sapovirus")
ictv_sapovirus_full<-subset(ictv_sapovirus,type=="full")
ictv_sapovirus_partial<-subset(ictv_sapovirus,type=="partial")
ictv_sapelovirus<-subset(ictv,molecule=="Sapelovirus")
ictv_sapelovirus_full<-subset(ictv_sapelovirus,type=="full")
ictv_sapelovirus_partial<-subset(ictv_sapelovirus,type=="partial")
ictv_teschovirus<-subset(ictv,molecule=="Teschovirus")

ictv_batpicornavirus_pep<-subset(ictv_pep,molecule=="Bat picornavirus")
ictv_hepatovirus_pep<-subset(ictv_pep, molecule=="Hepatovirus")
ictv_kobuvirus_pep<-subset(ictv_pep,molecule=="Kobuvirus")
ictv_kunsagivirus_pep<-subset(ictv_pep,molecule=="Kunsagivirus")
ictv_mischivirus_pep<-subset(ictv_pep,molecule=="Mischivirus")
ictv_sapovirus_pep<-subset(ictv_pep, molecule=="Sapovirus")
ictv_sapovirus_full_pep<-subset(ictv_sapovirus_pep,type=="full")
ictv_sapovirus_partial_pep<-subset(ictv_sapovirus_pep,type=="partial")
ictv_sapelovirus_pep<-subset(ictv_pep,molecule=="Sapelovirus")
ictv_sapelovirus_pep_full<-subset(ictv_sapelovirus_pep,type=="full")
ictv_sapelovirus_pep_partial<-subset(ictv_sapelovirus_pep,type=="partial")
ictv_teschovirus_pep<-subset(ictv_pep,molecule=="Teschovirus")

ictv_batpicornavirus_feat<-subset(ictv_feat,molecule=="Bat picornavirus")
ictv_hepatovirus_feat<-subset(ictv_feat,molecule=="Hepatovirus")
ictv_kobuvirus_feat<-subset(ictv_feat,molecule=="Kobuvirus")
ictv_kunsagivirus_feat<-subset(ictv_feat,molecule=="Kunsagivirus")
ictv_mischivirus_feat<-subset(ictv_feat,molecule=="Mischivirus")
ictv_sapovirus_feat<-subset(ictv_feat,molecule=="Sapovirus")
ictv_sapovirus_full_feat<-subset(ictv_sapovirus_feat,type=="full")
ictv_sapovirus_partial_feat<-subset(ictv_sapovirus_feat,type=="partial")
ictv_sapelovirus_feat<-subset(ictv_feat,molecule=="Sapelovirus")
ictv_sapelovirus_full_feat<-subset(ictv_sapelovirus_feat,type=="full")
ictv_sapelovirus_feat_partial<-subset(ictv_sapelovirus_feat,type=="partial")
ictv_teschovirus_feat<-subset(ictv_feat,molecule=="Teschovirus")

#plot ictv and blast plots
ictv_batpicorna<-ggplot(ictv_batpicornavirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_batpicornavirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                  xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_batpicornavirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,8200),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_batpicorna

ictv_hepato<-ggplot(ictv_hepatovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_hepatovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_hepatovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,6300),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_hepato


ictv_kobu<-ggplot(ictv_kobuvirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_kobuvirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_kobuvirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_kobuvirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                            xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_kobuvirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,9377),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_kobu


ictv_kun<-ggplot(ictv_kunsagivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_kunsagivirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_kunsagivirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_kunsagivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                               xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_kunsagivirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,7950),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_kun


ictv_mischi<-ggplot(ictv_mischivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_mischivirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_mischivirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_mischivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_mischivirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,9000),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_mischi


ictv_sapelo_full<-ggplot(ictv_sapelovirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_sapelovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_sapelovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_sapelovirus_pep_full, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                   xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_sapelovirus_full_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,8350),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapelo_full


ictv_tescho<-ggplot(ictv_teschovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_teschovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_teschovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_teschovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_teschovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,7400),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_tescho


ictv_sapo_full<-ggplot(ictv_sapovirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_sapovirus_full_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_sapovirus_full_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_sapovirus_full_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_sapovirus_full_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1000,7700),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapo_full



##Now get the plots from bootscan

#First all bat picornavirales over 3kb
setwd("~/Desktop/developer/mada-bat-picornavirus/bootscan/output_bat/all")

colzpalette<-c("#8ECAE6","#219EBC","#023047","#FFB703","#FB8500","#E48B97","#B52B09","#A60067","#987B6F","#8FD694")

#Bat picornavirus 
batpicorna_bat_boot <- read.csv(file = "batpicornavirus_boot.csv", header = T, stringsAsFactors = F)
head(batpicorna_bat_boot)

#move to long
long.sim_nt <- melt(batpicorna_bat_boot, id.vars = c("pointer"), measure.vars = c("OQ818328","HQ595345","KJ641687","HQ595344", "KJ641690","KJ641699","MF352419","MF352423"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818328"] <- "Rousettus madagascariensis picornavirus OQ818328"
long.sim_nt$accession[long.sim_nt$accession == "HQ595345"] <- "Hipposideros armiger picornavirus HQ595345"
long.sim_nt$accession[long.sim_nt$accession == "KJ641687"] <- "Miniopterus fuliginosus picornavirus KJ641687"
long.sim_nt$accession[long.sim_nt$accession == "HQ595344"] <- "Rhinolophus sinicus picornavirus HQ595344"
long.sim_nt$accession[long.sim_nt$accession == "KJ641690"] <- "Miniopterus fuliginosus picornavirus KJ641690"
long.sim_nt$accession[long.sim_nt$accession == "KJ641699"] <- "Miniopterus fuliginosus picornavirus KJ641699"
long.sim_nt$accession[long.sim_nt$accession == "MF352419"] <- "Miniopterus schreibersii picornavirus MF352419"
long.sim_nt$accession[long.sim_nt$accession == "MF352423"] <- "Miniopterus schreibersii picornavirus MF352423"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rousettus madagascariensis picornavirus OQ818328", "Hipposideros armiger picornavirus HQ595345",
                                                                  "Miniopterus fuliginosus picornavirus KJ641687","Rhinolophus sinicus picornavirus HQ595344",
                                                                  "Miniopterus fuliginosus picornavirus KJ641690","Miniopterus fuliginosus picornavirus KJ641699",
                                                                  "Miniopterus schreibersii picornavirus MF352419",
                                                                  "Miniopterus schreibersii picornavirus MF352423"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("Reference: ",italic("Rousettus madagascariensis picornavirus "), "OQ818325"))

batpicorna_bat_boot <- ggplot(long.sim_nt) + 
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
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
  #                    labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_bat_boot

#put gene map with PySimPlot
batpicorna_boot<-batpicorna_bat_boot/bat_batpicorna_all+plot_layout(nrow=2,  heights = c(2, 0.30))
batpicorna_boot

batpicorna_boot<-as.ggplot(batpicorna_boot)
batpicorna_boot



#Hepatovirus
hepatovirus_bat_all_boot <- read.csv(file = "hepatovitus_boot.csv", header = T, stringsAsFactors = F) #animo acid
head(hepatovirus_bat_all_boot)

#move to long
long.sim_nt <- melt(hepatovirus_bat_all_boot, id.vars = c("pointer"), measure.vars = c("KT452714","KX420952","MG559674","KT452730",
                                                                                     "KT452729","OM302498","KT452742"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "KT452714"] <- "Eidolon helvum hepatovirus KT452714"
long.sim_nt$accession[long.sim_nt$accession == "KX420952"] <- "Rhinopoma hardwickii hepatovirus KX420952"
long.sim_nt$accession[long.sim_nt$accession == "MG559674"] <- "Hipposideros armiger hepatovirus MG559674"
long.sim_nt$accession[long.sim_nt$accession == "KT452729"] <- "Rhinolophus landeri hepatovirus KT452729"
long.sim_nt$accession[long.sim_nt$accession == "KT452730"] <- "Coleura afra hepatovirus KT452730"
long.sim_nt$accession[long.sim_nt$accession == "OM302498"] <- "Eptesicus fuscus hepatovirus OM302498"
long.sim_nt$accession[long.sim_nt$accession == "KT452742"] <- "Miniopterus cf. manavi hepatovirus KT452742"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon helvum hepatovirus KT452714",
                                                                  "Rhinopoma hardwickii hepatovirus KX420952",
                                                                  "Hipposideros armiger hepatovirus MG559674",
                                                                  "Rhinolophus landeri hepatovirus KT452729",
                                                                  "Coleura afra hepatovirus KT452730",
                                                                  "Eptesicus fuscus hepatovirus OM302498",
                                                                  "Miniopterus cf. manavi hepatovirus KT452742"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum hepatovirus "), "OQ818337"))

hepatovirus_bat_all_boot <- ggplot(long.sim_nt) + 
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
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_bat_all_boot


#put gene map with PySimPlot
hep_bat_all_boot<-hepatovirus_bat_all_boot/bat_hepato_all+plot_layout(nrow=2,  heights = c(2, 0.30))
hep_bat_all_boot

hep_bat_all_boot<-as.ggplot(hep_bat_all_boot)
hep_bat_all_boot



#mischivirus
mischivirus_bat_all_boot <- read.csv(file = "mischivirus_boot.csv", header = T, stringsAsFactors = F) #Nucleotide
head(mischivirus_bat_all_boot)

#move to long
long.sim_nt <- melt(mischivirus_bat_all_boot, id.vars = c("pointer"), measure.vars = c("JQ814851","KP054273",
                                                                                "KP054276","KP054277",
                                                                                "KP054278",
                                                                                "MG888045","KP100644"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "JQ814851"] <- "Miniopterus schreibersii mischivirus JQ814851"
long.sim_nt$accession[long.sim_nt$accession == "KP054273"] <- "Miniopterus schreibersii mischivirus KP054273"
long.sim_nt$accession[long.sim_nt$accession == "KP054276"] <- "Miniopterus schreibersii shanbavirus KP054276"
long.sim_nt$accession[long.sim_nt$accession == "KP054277"] <- "Miniopterus schreibersii shanbavirus KP054277"
long.sim_nt$accession[long.sim_nt$accession == "KP054278"] <- "Miniopterus schreibersii shanbavirus KP054278"
long.sim_nt$accession[long.sim_nt$accession == "MG888045"] <- "Miniopterus sp. mischivirus MG888045"
long.sim_nt$accession[long.sim_nt$accession == "KP100644"] <- "Hipposideros gigas mischivirus KP100644"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Miniopterus schreibersii mischivirus JQ814851",
                                                                  "Miniopterus schreibersii mischivirus KP054273",
                                                                  "Miniopterus schreibersii shanbavirus KP054276",
                                                                  "Miniopterus schreibersii shanbavirus KP054277",
                                                                  "Miniopterus schreibersii shanbavirus KP054278",
                                                                  "Miniopterus sp. mischivirus MG888045",
                                                                  "Hipposideros gigas mischivirus KP100644"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Pteropus rufus mischivirus "), "OQ818316"))

mischivirus_bat_all_boot <- ggplot(long.sim_nt) + 
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
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_bat_all_boot

#put gene map with PySimPlot
mischi_bat_all_boot<-mischivirus_bat_all_boot/bat_mischi_all+plot_layout(nrow=2,  heights = c(2, 0.30))
mischi_bat_all_boot

mischi_bat_all_boot<-as.ggplot(mischi_bat_all_boot)
mischi_bat_all_boot


#Sapelovirus
sapelovirus_bat_all_boot <- read.csv(file = "sapelovirus_boot.csv", header = T, stringsAsFactors = F) #Nucleotide
head(sapelovirus_bat_all_boot)

#move to long
long.sim_nt <- melt(sapelovirus_bat_all_boot, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","KJ641689","KJ641696",
                                                                                       "MF352431","KX644938"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_nt$accession[long.sim_nt$accession == "KJ641696"] <- "Vespertilio superans shanbavirus KJ641696"
long.sim_nt$accession[long.sim_nt$accession == "KJ641689"] <- "Myotis altarium shanbavirus KJ641689"
long.sim_nt$accession[long.sim_nt$accession == "MF352431"] <- "Myotis ricketti sapelovirus MF352431"
long.sim_nt$accession[long.sim_nt$accession == "KX644938"] <- "Eidolon helvum sapelovirus KX644938"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321",
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Vespertilio superans shanbavirus KJ641696",
                                                                  "Myotis altarium shanbavirus KJ641689",
                                                                  "Myotis ricketti sapelovirus MF352431",
                                                                  "Eidolon helvum sapelovirus KX644938"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapelovirus "), "OQ818320"))

sapelovirus_bat_all_boot <- ggplot(long.sim_nt) + 
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
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_bat_all_boot

#put gene map with PySimPlot
sapelo_bat_all_boot<-sapelovirus_bat_all_boot/bat_sapelo_all+plot_layout(nrow=2,  heights = c(2, 0.30))
sapelo_bat_all_boot

sapelo_bat_all_boot<-as.ggplot(sapelo_bat_all_boot)
sapelo_bat_all_boot




#Sapovirus
sapovirus_all_bat_boot <- read.csv(file = "sapovirus_boot.csv", header = T, stringsAsFactors = F) #Nucleotide
head(sapovirus_all_bat_boot)

#move to long
long.sim_nt <- melt(sapovirus_all_bat_boot, id.vars = c("pointer"), measure.vars = c("KX759618","KX759620","KX759621","KX759622",
                                                                                    "KX759619", "KX759623","OP963623"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "KX759618"] <- "Eidolon helvum sapovirus KX759618"
long.sim_nt$accession[long.sim_nt$accession == "KX759620"] <- "Eidolon helvum sapovirus KX759620"
long.sim_nt$accession[long.sim_nt$accession == "KX759621"] <- "Eidolon helvum sapovirus KX759621"
long.sim_nt$accession[long.sim_nt$accession == "KX759622"] <- "Eidolon helvum sapovirus KX759622"
long.sim_nt$accession[long.sim_nt$accession == "KX759619"] <- "Eidolon helvum sapovirus KX759619"
long.sim_nt$accession[long.sim_nt$accession == "KX759623"] <- "Eidolon helvum sapovirus KX759623"
long.sim_nt$accession[long.sim_nt$accession == "OP963623"] <- "Rousettus leschenaultii OP963623"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon helvum sapovirus KX759618",
                                                                  "Eidolon helvum sapovirus KX759620",
                                                                  "Eidolon helvum sapovirus KX759621",
                                                                  "Eidolon helvum sapovirus KX759622",
                                                                  "Eidolon helvum sapovirus KX759619",
                                                                  "Eidolon helvum sapovirus KX759623", 
                                                                  "Rousettus leschenaultii OP963623"))
#and plot
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapovirus "), "OQ818319"))

sapovirus_bat_boot <- ggplot(long.sim_nt) + 
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
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_bat_boot

#put gene map with PySimPlot
sapo_bat_boot<-sapovirus_bat_boot/bat_sapo_all+plot_layout(nrow=2,  heights = c(2, 0.30))
sapo_bat_boot

sapo_bat_boot<-as.ggplot(sapo_bat_boot)
sapo_bat_boot



#Teschovirus
teschovirus_bat_boot <- read.csv(file = "teschovirus_boot.csv", header = T, stringsAsFactors = F) #Nucleotide
head(teschovirus_bat_boot)

#move to long
long.sim_nt <- melt(teschovirus_bat_boot, id.vars = c("pointer"), measure.vars = c("OQ818323","OQ818324","KX420938"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818323"] <- "Rousettus madagascariensis teschovirus OQ818323"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324"] <- "Rousettus madagascariensis teschovirus OQ818324"
long.sim_nt$accession[long.sim_nt$accession == "KX420938"] <- "Eidolon helvum teschovirus KX420938"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rousettus madagascariensis teschovirus OQ818323",
                                                                  "Rousettus madagascariensis teschovirus OQ818324",
                                                                  "Eidolon helvum teschovirus KX420938"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum teschovirus "), "OQ818318"))

teschovirus_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
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
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

teschovirus_bat_boot

#put gene map with PySimPlot
tescho_boot<-teschovirus_bat_boot/bat_tescho_all+plot_layout(nrow=2,  heights = c(2, 0.30))
tescho_boot

tescho_boot<-as.ggplot(tescho_boot)
tescho_boot










##Now plot the bootscans using only full picornavirales as reference
setwd("~/Desktop/developer/mada-bat-picornavirus/bootscan/output_bat/full")


#Hepatovirus
hepatovirus_bat_full_boot <- read.csv(file = "hepatovirus_boot.csv", header = T, stringsAsFactors = F) #animo acid
head(hepatovirus_bat_full_boot)

#move to long
long.sim_nt <- melt(hepatovirus_bat_full_boot, id.vars = c("pointer"), measure.vars = c("NC_028366","NC_038316","NC_038313"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_028366"] <- "Eidolon helvum hepatovirus NC_028366"
long.sim_nt$accession[long.sim_nt$accession == "NC_038316"] <- "Coleura afra hepatovirus NC_038316"
long.sim_nt$accession[long.sim_nt$accession == "NC_038313"] <- "Miniopterus cf. manavi hepatovirus NC_038313"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon helvum hepatovirus NC_028366",
                                                                  "Coleura afra hepatovirus NC_038316",
                                                                  "Miniopterus cf. manavi hepatovirus NC_038313"))

long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum hepatovirus "), "OQ818337"))

hepatovirus_bat_full_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
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
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_bat_full_boot


#put gene map with PySimPlot
hep_bat_full_boot<-hepatovirus_bat_full_boot/bat_hepato_full+plot_layout(nrow=2,  heights = c(2, 0.30))
hep_bat_full_boot

hep_bat_full_boot<-as.ggplot(hep_bat_full_boot)
hep_bat_full_boot




#kunsagivirus
kunsagivirus_bat_full_boot <- read.csv(file = "kunsagivirus_boot.csv", header = T, stringsAsFactors = F) #Nucleotide
head(kunsagivirus_bat_full_boot)

#move to long
long.sim_nt <- melt(kunsagivirus_bat_full_boot, id.vars = c("pointer"), measure.vars = c("MN602325","NC_043071",
                                                                                        "NC_033818"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "MN602325"] <- "Pipistrellus pipistrellus shanbavirus MN602325"
long.sim_nt$accession[long.sim_nt$accession == "NC_043071"] <- "Myotis ricketti shanbavirus NC_043071"
long.sim_nt$accession[long.sim_nt$accession == "NC_033818"] <- "Eidolon helvum kunsagivirus NC_033818"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Pipistrellus pipistrellus shanbavirus MN602325",
                                                                  "Myotis ricketti shanbavirus NC_043071",
                                                                  "Eidolon helvum kunsagivirus NC_033818"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum kunsagivirus "), "OQ818317"))

kunsagivirus_bat_full_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
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
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagivirus_bat_full_boot

#put gene map with PySimPlot
kun_bat_full_boot<-kunsagivirus_bat_full_boot/bat_kun_full+plot_layout(nrow=2,  heights = c(2, 0.30))
kun_bat_full_boot

kun_bat_full_boot<-as.ggplot(kun_bat_full_boot)
kun_bat_full_boot



#mischivirus
mischivirus_bat_full_boot <- read.csv(file = "mischivirus_boot.csv", header = T, stringsAsFactors = F) #Nucleotide
head(mischivirus_bat_full_boot)

#move to long
long.sim_nt <- melt(mischivirus_bat_full_boot, id.vars = c("pointer"), measure.vars = c("JQ814851","NC_043072","KP100644"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "JQ814851"] <- "Miniopterus schreibersii mischivirus JQ814851"
long.sim_nt$accession[long.sim_nt$accession == "NC_043072"] <-"Miniopterus schreibersii mischivirus NC_043072"
long.sim_nt$accession[long.sim_nt$accession == "KP100644"] <- "Hipposideros gigas mischivirus KP100644"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Miniopterus schreibersii mischivirus JQ814851",
                                                                  "Miniopterus schreibersii mischivirus NC_043072",
                                                                  "Hipposideros gigas mischivirus KP100644"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Pteropus rufus mischivirus "), "OQ818316"))

mischivirus_bat_full_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
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
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_bat_full_boot

#put gene map with PySimPlot
mischi_bat_full_boot<-mischivirus_bat_full_boot/bat_mischi_full+plot_layout(nrow=2,  heights = c(2, 0.30))
mischi_bat_full_boot

mischi_bat_full_boot<-as.ggplot(mischi_bat_full_boot)
mischi_bat_full_boot


#Sapelovirus
sapelovirus_bat_full_boot <- read.csv(file = "sapelovirus_boot.csv", header = T, stringsAsFactors = F) #Nucleotide
head(sapelovirus_bat_full_boot)

#move to long
long.sim_nt <- melt(sapelovirus_bat_full_boot, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","NC_033820",
                                                                                        "NC_076039","NC_076040"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820"] <- "Eidolon helvum sapelovirus NC_033820"
long.sim_nt$accession[long.sim_nt$accession == "NC_076039"] <- "Rhinolophus ferrumequinum shanbavirus NC_076039"
long.sim_nt$accession[long.sim_nt$accession == "NC_076040"] <- "Nyctalus velutinus shanbavirus NC_076040"



long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321",
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Eidolon helvum sapelovirus NC_033820",
                                                                  "Rhinolophus ferrumequinum shanbavirus NC_076039",
                                                                  "Nyctalus velutinus shanbavirus NC_076040"
                                                                  
))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapelovirus "), "OQ818320"))

sapelovirus_bat_full_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
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
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_bat_full_boot

#put gene map with PySimPlot
sapelo_bat_full_boot<-sapelovirus_bat_full_boot/bat_sapelo_full+plot_layout(nrow=2,  heights = c(2, 0.30))
sapelo_bat_full_boot

sapelo_bat_full_boot<-as.ggplot(sapelo_bat_full_boot)
sapelo_bat_full_boot





##Now put the whole figure together

#all bootscan together: all bat picornavirales
bootscan_batall<-plot_grid(batpicorna_boot,
                           hep_bat_all_boot,
                           mischi_bat_all_boot,
                           sapelo_bat_all_boot,
                           sapo_bat_boot,
                           tescho_boot,
                    ncol=3,
                    labels="AUTO",  label_size = 23, align = "hv", axis="b")
bootscan_batall


#all bootscan together: full bat picornavirales reference
bootscan_batfull<-plot_grid(hep_bat_full_boot,
                           kun_bat_full_boot,
                           mischi_bat_full_boot,
                           sapelo_bat_full_boot,
                           ncol=2,
                           labels="AUTO",  label_size = 23, align = "hv", axis="b")
bootscan_batfull





#in Fig 3
bootscan_bat_fig<-plot_grid(mischi_bat_all_boot,
                            sapelo_bat_all_boot,
                            hep_bat_all_boot,
                        nrow=1,
                        labels="AUTO",  label_size = 23, align = "hv", axis="b")
bootscan_bat_fig


#excluded from Fig 3 
bootscan_bat_supp<-plot_grid(batpicorna_boot,
                             sapo_bat_boot,
                             kun_bat_full_boot,
                             mischi_bat_full_boot,
                             sapelo_bat_full_boot,
                             hep_bat_full_boot,
                         ncol=3,
                         labels="AUTO",  label_size = 23, align = "hv", axis="b")
bootscan_bat_supp

