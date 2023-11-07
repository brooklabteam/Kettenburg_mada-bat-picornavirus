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


#This is to make a figure of the representative PySimPlots with their matching genome plots

##First make the genome plots
homewd="/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/genome_annotation_and_characterization/"
setwd("~/Desktop/developer/mada-bat-picornavirus/genome_annotation_and_characterization/genus_gene_maps")

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
# colz=c("5'UTR"="gold", "L"="royalblue","VP4"="paleturquoise3", "VP2"="skyblue1", "VP0"="royalblue4", "VP3"="steelblue1",
#        "VP1"="cadetblue1", "VP1/2A"="cadetblue1", "2A"="orange1", "2B"="sienna2", "2C"="darkorange1", "3A"="palevioletred1", "3B"="plum",
#        "3C"="rosybrown1", "3D"="pink2", "Helicase"="darkseagreen1","NS4"="darkolivegreen1","Vpg"="seagreen1","Pro-Pol"="palegreen2",
#        "Polyprotein"="azure3","Putative polyprotein"="mediumorchid1", "Non-structural polyprotein"="mediumorchid4",
#        "Putative minor structural protein" ="slateblue3", "Structural polyprotein"="lightgoldenrod1",
#        "Hypothetical protein"="darkslategrey", "Similar to structural polyprotein"="mediumpurple1", "Similar to putative polyprotein"="mediumpurple3",
#        "Similar to polyprotein"="mediumpurple4","3'UTR"="yellow")

colz=c("5'UTR"="grey", "L"="white","VP4"="white", "VP2"="white", "VP0"="white", "VP3"="white",
       "VP1"="white", "VP1/2A"="white", "2A"="white", "2B"="white", "2C"="white", "3A"="white", "3B"="white",
       "3C"="white", "3D"="white", "Helicase"="white","NS4"="white","Vpg"="white","Pro-Pol"="white",
       "Polyprotein"="black","Putative polyprotein"="black", "Non-structural polyprotein"="black",
       "Putative minor structural protein" ="white", "Structural polyprotein"="black",
       "Hypothetical protein"="white", "Similar to structural polyprotein"="black", "Similar to putative polyprotein"="black",
       "Similar to polyprotein"="black","3'UTR"="grey")


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

#gene maps
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
                                                                  xsubmin=from, xsubmax=to), color="black",
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
                                                              xsubmin=from, xsubmax=to), color="black",
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
                                                            xsubmin=from, xsubmax=to), color="black",
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
                                                               xsubmin=from, xsubmax=to), color="black",
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
                                                              xsubmin=from, xsubmax=to), color="black",
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
                                                                   xsubmin=from, xsubmax=to), color="black",
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
                                                              xsubmin=from, xsubmax=to), color="black",
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
                                                                 xsubmin=from, xsubmax=to), color="black",
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



##Now get the plots for the PySimPlot, just the ICTV_BLAST ones, there will be a separate PDF file with the table of African bat picorna similarities
setwd("~/Desktop/developer/mada-bat-picornavirus/PySimPlot/ICTV_BLAST_pysimplot")

colzpalette<-c("coral1","cyan3","hotpink1","gold","orange3","dodgerblue2","firebrick1","darkorchid3","skyblue1")


#Bat picornavirus
bat_picorna_ictv_aa_full <- read.csv(file = "bat_picorna_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F)
head(bat_picorna_ictv_aa_full)

#move to long
long.sim_aa <- melt(bat_picorna_ictv_aa_full, id.vars = c("pointer"), measure.vars = c("OQ818328","HQ595340","HQ595342","HQ595344"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "OQ818328"] <- "R. madagascariensis roupivirus MIZ240"
long.sim_aa$accession[long.sim_aa$accession == "HQ595340"] <- "Bat picornavirus 1: HQ595340"
long.sim_aa$accession[long.sim_aa$accession == "HQ595342"] <- "Bat picornavirus 2: HQ595342"
long.sim_aa$accession[long.sim_aa$accession == "HQ595344"] <- "Bat picornavirus 3: HQ595344"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("R. madagascariensis roupivirus MIZ240", "Bat picornavirus 1: HQ595340",
                                                                  "Bat picornavirus 2: HQ595342","Bat picornavirus 3: HQ595344"))
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

title<-expression(paste("Reference: R. madagascariensis roupivirus MIZ214"))

## Nucleotide
batpicorna_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
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
  ggtitle(title)+
  scale_color_manual(values=colzpalette) + 
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_ictv_aa

#put gene map with PySimPlot
bat_picorna_aa<-batpicorna_ictv_aa/ictv_batpicorna+plot_layout(nrow=2,  heights = c(1, 0.3))
bat_picorna_aa

bat_picorna_aa<-as.ggplot(bat_picorna_aa)
bat_picorna_aa



#Hepatovirus
hepato_ictv_aa_partial <- read.csv(file = "hepato_ictv_aa_partial_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(hepato_ictv_aa_partial)

#move to long
long.sim_aa <- melt(hepato_ictv_aa_partial, id.vars = c("pointer"), measure.vars = c("HPA","KR703607","KT452637","KT452658",
                                                                                     "KT452685","KT452714","KT452730", "KT452735",
                                                                                     "KT452742"))
unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "HPA"] <- "Hepatovirus A HPA"
long.sim_aa$accession[long.sim_aa$accession == "KR703607"] <- "Hepatovirus B: KR703607"
long.sim_aa$accession[long.sim_aa$accession == "KT452637"] <- "Hepatovirus D: KT452637"
long.sim_aa$accession[long.sim_aa$accession == "KT452658"] <- "Hepatovirus I: KT452658"
long.sim_aa$accession[long.sim_aa$accession == "KT452685"] <- "Hepatovirus F: KT452685"
long.sim_aa$accession[long.sim_aa$accession == "KT452714"] <- "Hepatovirus H: KT452714"
long.sim_aa$accession[long.sim_aa$accession == "KT452730"] <- "Hepatovirus G: KT452730"
long.sim_aa$accession[long.sim_aa$accession == "KT452735"] <- "Hepatovirus E: KT452735"
long.sim_aa$accession[long.sim_aa$accession == "KT452742"] <- "Hepatovirus C: KT452742"


long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Hepatovirus A HPA", "Hepatovirus B: KR703607",
                                                                  "Hepatovirus C: KT452742","Hepatovirus D: KT452637",
                                                                  "Hepatovirus E: KT452735","Hepatovirus F: KT452685",
                                                                  "Hepatovirus G: KT452730", "Hepatovirus H: KT452714",
                                                                  "Hepatovirus I: KT452658"))
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## Amino acid
title<-expression(paste("Reference: E. dupreanum hepatovirus KEL148"))

hepatovirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        #legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  #scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055),
                     labels = c(0,2000, 4000,6000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_ictv_aa


#put gene map with PySimPlot
hep_ictv_aa<-hepatovirus_ictv_aa/ictv_hepato+plot_layout(nrow=2,  heights = c(1, 0.3))
hep_ictv_aa

hep_ictv_aa<-as.ggplot(hep_ictv_aa)
hep_ictv_aa


#kunsagivirus
kunsagivirus_aa_ictv <- read.csv(file = "kun_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #Amino acid
head(kunsagivirus_aa_ictv)

#move to long
long.sim_aa <- melt(kunsagivirus_aa_ictv, id.vars = c("pointer"), measure.vars = c("KC935379","KX644936",
                                                                                   "KY670597"))
unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "KC935379"] <- "Kunsagivirus A: KC935379"
long.sim_aa$accession[long.sim_aa$accession == "KX644936"] <- "Kunsagivirus B: KX644936"
long.sim_aa$accession[long.sim_aa$accession == "KY670597"] <- "Kunsagivirus C: KY670597"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Kunsagivirus A: KC935379", "Kunsagivirus B: KX644936",
                                                                  "Kunsagivirus C: KY670597"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## Amino acid
title<-expression(paste("Reference: E. dupreanum kunsagivirus KEL148"))

kunsagivirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        # legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  #scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagivirus_ictv_aa

#put gene map with PySimPlot
kun_ictv_aa<-kunsagivirus_ictv_aa/ictv_kun+plot_layout(nrow=2,  heights = c(1, 0.3))
kun_ictv_aa

kun_ictv_aa<-as.ggplot(kun_ictv_aa)
kun_ictv_aa


#Mischivirus
mischivirus_aa_ictv <- read.csv(file = "mischi_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #Amino acid
head(mischivirus_aa_ictv)

#move to long
long.sim_aa <- melt(mischivirus_aa_ictv, id.vars = c("pointer"), measure.vars = c("JQ814851","KP054273",
                                                                                  "KP100644","KY512802",
                                                                                  "MF352410"))
unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "JQ814851"] <- "Mischivirus A: JQ814851"
long.sim_aa$accession[long.sim_aa$accession == "KP054273"] <- "Mischivirus B: KP054273"
long.sim_aa$accession[long.sim_aa$accession == "KP100644"] <- "Mischivirus C: KP100644"
long.sim_aa$accession[long.sim_aa$accession == "KY512802"] <- "Mischivirus D: KY512802"
long.sim_aa$accession[long.sim_aa$accession == "MF352410"] <- "Mischivirus E: MF352410"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Mischivirus A: JQ814851", "Mischivirus B: KP054273",
                                                                  "Mischivirus C: KP100644", "Mischivirus D: KY512802",
                                                                  "Mischivirus E: MF352410"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## Amino acid
title<-expression(paste("Reference: P. rufus mischivirus AMB150"))

mischivirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        #legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  #scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_ictv_aa

#put gene map with PySimPlot
mischi_ictv_aa<-mischivirus_ictv_aa/ictv_mischi+plot_layout(nrow=2,  heights = c(1, 0.3))
mischi_ictv_aa

mischi_ictv_aa<-as.ggplot(mischi_ictv_aa)
mischi_ictv_aa



#Sapelovirus full
sapelovirus_full_aa_ictv <- read.csv(file = "sapelo_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapelovirus_full_aa_ictv)

#move to long
long.sim_aa <- melt(sapelovirus_full_aa_ictv, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","AF406813","AY064708",
                                                                                       "NC_033820"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "OQ818321"] <- "E. dupreanum sapelovirus KEL272"
long.sim_aa$accession[long.sim_aa$accession == "OQ818329"] <- "R. madagascariensis sapelovirus MIZ243"
long.sim_aa$accession[long.sim_aa$accession == "AF406813"] <- "Sapelovirus A: AF406813"
long.sim_aa$accession[long.sim_aa$accession == "AY064708"] <- "Sapelovirus B: AY064708"
long.sim_aa$accession[long.sim_aa$accession == "NC_033820"] <- "Eidolon helvum sapelovirus: NC_033820"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("E. dupreanum sapelovirus KEL272",
                                                                  "R. madagascariensis sapelovirus MIZ243",
                                                                  "Sapelovirus A: AF406813","Sapelovirus B: AY064708",
                                                                  "Eidolon helvum sapelovirus: NC_033820"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## Amino acid
title<-expression(paste("Reference: E. dupreanum sapelovirus KEL233"))

sapelovirus_full_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        #legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  #scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_full_ictv_aa

#put gene map with PySimPlot
sapelo_full_ictv_aa<-sapelovirus_full_ictv_aa/ictv_sapelo_full+plot_layout(nrow=2,  heights = c(1, 0.3))
sapelo_full_ictv_aa

sapelo_full_ictv_aa<-as.ggplot(sapelo_full_ictv_aa)
sapelo_full_ictv_aa


#Teschovirus
teschovirus_aa_ictv <- read.csv(file = "tescho_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_aa_ictv)

#move to long
long.sim_aa <- melt(teschovirus_aa_ictv, id.vars = c("pointer"), measure.vars = c("OQ818323","OQ818324","LC386158","MG875515","MT295502"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "OQ818323"] <- "R. madagascariensis teschovirus MIZ190"
long.sim_aa$accession[long.sim_aa$accession == "OQ818324"] <- "R. madagascariensis teschovirus MIZ205"
long.sim_aa$accession[long.sim_aa$accession == "LC386158"] <- "Teschovirus A: LC386158"
long.sim_aa$accession[long.sim_aa$accession == "MG875515"] <- "Teschovirus B: MG875515"
long.sim_aa$accession[long.sim_aa$accession == "MT295502"] <- "Teschovirus A6: MT295502"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("R. madagascariensis teschovirus MIZ190",
                                                                  "R. madagascariensis teschovirus MIZ205",
                                                                  "Teschovirus A: LC386158",
                                                                  "Teschovirus B: MG875515",
                                                                  "Teschovirus A6: MT295502"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## Amino acid
title<-expression(paste("Reference: E. dupreanum teschovirus KEL164"))

teschovirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        #legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  #scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

teschovirus_ictv_aa

#put gene map with PySimPlot
tescho_ictv_aa<-teschovirus_ictv_aa/ictv_tescho+plot_layout(nrow=2,  heights = c(1, 0.3))
tescho_ictv_aa

tescho_ictv_aa<-as.ggplot(tescho_ictv_aa)
tescho_ictv_aa



#Sapovirus full
sapovirus_full_nt_ictv <- read.csv(file = "sapo_ictv_nt_full_alignment_trimmed.csv", header = T, stringsAsFactors = F) #nucleotide
head(sapovirus_full_nt_ictv)

#move to long
long.sim_nt <- melt(sapovirus_full_nt_ictv, id.vars = c("pointer"), measure.vars = c("HM002617","KX759623"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "HM002617"] <- "Sapporo virus: HM002617"
long.sim_nt$accession[long.sim_nt$accession == "KX759623"] <- "Bat sapovirus: KX759623"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Sapporo virus: HM002617", "Bat sapovirus: KX759623"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste("E. dupreanum sapovirus KEL166"))

sapovirus_full_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face="italic", size = 8),
        #legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
  #                    labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_full_ictv_nt


sapovirus_full_aa_ictv <- read.csv(file = "sapo_ictv_aa_full_alignment_trimmed.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapovirus_full_aa_ictv)

#move to long
long.sim_aa <- melt(sapovirus_full_aa_ictv, id.vars = c("pointer"), measure.vars = c("HM002617","KX759623"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "HM002617"] <- "Sapporo virus: HM002617"
long.sim_aa$accession[long.sim_aa$accession == "KX759623"] <- "Bat sapovirus: KX759623"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Sapporo virus: HM002617", "Bat sapovirus: KX759623"))

#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## Amino acid
sapovirus_full_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="none", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        #legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_full_ictv_aa

#put gene map with PySimPlot
sapo_full_ictv<-sapovirus_full_ictv_nt/sapovirus_full_ictv_aa/ictv_sapo_full+plot_layout(nrow=3,  heights = c(0.5,0.5, 0.3))
sapo_full_ictv

sapo_full_ictv<-as.ggplot(sapo_full_ictv)
sapo_full_ictv

#export 6x8 landscape PDF for sapovirus only figure 3


##Now put the whole figure together
fig2<-plot_grid(mischi_ictv_aa,bat_picorna_aa,hep_ictv_aa,kun_ictv_aa,
                sapelo_full_ictv_aa,tescho_ictv_aa,
                ncol=3,
                labels="AUTO", label_size = 23, align = "hv", axis="b")
fig2


#export 13x23 landscape PDF

