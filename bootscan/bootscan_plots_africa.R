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
ictv <- read.csv("ictv_blast_genes_trimmed.csv", header = T, stringsAsFactors = F)
ictv$gene<-factor(ictv$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                        "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                        "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                        "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                        "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the peptide files
ictv_pep <- read.csv("ictv_blast_peptides_trimmed.csv", header = T, stringsAsFactors = F)
ictv_pep$gene<-factor(ictv_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                                "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                                "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the feature file in case its needed
ictv_feat <- read.csv("ictv_blast_features_trimmed.csv", header = T, stringsAsFactors = F)
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
  scale_x_continuous(limits=c(4540,5670),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,6380),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,8300),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,9350),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,8150),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,8950),expand=c(0,0))+
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
  scale_x_continuous(limits=c(2700,4350),expand=c(0,0))+
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
setwd("~/Desktop/developer/mada-bat-picornavirus/bootscan/output_africa")

colzpalette<-c("#8ECAE6","#219EBC","#023047","#FFB703","#FB8500","#E48B97","#B52B09","#A60067","#987B6F","#8FD694")

#Bat picornavirus 
africa_batpicorna_bootscan <- read.csv(file = "africa_batpicorna_bootscan.csv", header = T, stringsAsFactors = F)
head(africa_batpicorna_bootscan)

#move to long
long.sim_nt <- melt(africa_batpicorna_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818325","OQ818346","JX437642","KF040078", "KF040079"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818325"] <- "Rousettus madagascariensis picornavirus OQ818325"
long.sim_nt$accession[long.sim_nt$accession == "OQ818346"] <- "Rousettus madagascariensis picornavirus OQ818346"
long.sim_nt$accession[long.sim_nt$accession == "JX437642"] <- "Homo sapiens enterovirus JX437642"
long.sim_nt$accession[long.sim_nt$accession == "KF040078"] <- "Pan troglodytes enterovirus KF040078"
long.sim_nt$accession[long.sim_nt$accession == "KF040079"] <- "Pan troglodytes enterovirus KF040079"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rousettus madagascariensis picornavirus OQ818325", "Rousettus madagascariensis picornavirus OQ818346",
                                                                  "Homo sapiens enterovirus JX437642","Pan troglodytes enterovirus KF040078",
                                                                  "Pan troglodytes enterovirus KF040079"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("Reference: ",italic("Rousettus madagascariensis picornavirus "), "OQ818328"))

batpicorna_africa_boot <- ggplot(long.sim_nt) + 
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

batpicorna_africa_boot

#put gene map with PySimPlot
batpicorna_boot<-batpicorna_africa_boot/ictv_batpicorna+plot_layout(nrow=2,  heights = c(2, 0.30))
batpicorna_boot

batpicorna_boot<-as.ggplot(batpicorna_boot)
batpicorna_boot



#Hepatovirus
africa_hepato_bootscan <- read.csv(file = "africa_hepato_bootscan.csv", header = T, stringsAsFactors = F) #animo acid
head(africa_hepato_bootscan)

#move to long
long.sim_nt <- melt(africa_hepato_bootscan, id.vars = c("pointer"), measure.vars = c("AY644670","KT819575","KT452729","NC_038316",
                                                                                     "NC_028366","NC_038313","NC_038314"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "AY644670"] <- "Homo sapiens hepatovirus AY644670"
long.sim_nt$accession[long.sim_nt$accession == "KT819575"] <- "Papio anubis hepatovirus KT819575"
long.sim_nt$accession[long.sim_nt$accession == "KT452729"] <- "Rhinolophus landeri hepatovirus KT452729"
long.sim_nt$accession[long.sim_nt$accession == "NC_038316"] <- "Coleura afra hepatovirus hepatovirus NC_038316"
long.sim_nt$accession[long.sim_nt$accession == "NC_028366"] <- "Homo sapiens hepatovirus NC_028366"
long.sim_nt$accession[long.sim_nt$accession == "NC_038313"] <- "Miniopterus sp. hepatovirus NC_038313"
long.sim_nt$accession[long.sim_nt$accession == "NC_038314"] <- "Lophuromys sikapusi hepatovirus NC_038314"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Homo sapiens hepatovirus AY644670",
                                                                  "Papio anubis hepatovirus KT819575",
                                                                  "Rhinolophus landeri hepatovirus KT452729",
                                                                  "Coleura afra hepatovirus hepatovirus NC_038316",
                                                                  "Homo sapiens hepatovirus NC_028366",
                                                                  "Miniopterus sp. hepatovirus NC_038313",
                                                                  "Lophuromys sikapusi hepatovirus NC_038314"))
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
hep_bat_all_boot<-hepatovirus_bat_all_boot/ictv_hepato+plot_layout(nrow=2,  heights = c(2, 0.30))
hep_bat_all_boot

hep_bat_all_boot<-as.ggplot(hep_bat_all_boot)
hep_bat_all_boot



#mischivirus
africa_mischi_bootscan <- read.csv(file = "africa_mischi_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_mischi_bootscan)

#move to long
long.sim_nt <- melt(africa_mischi_bootscan, id.vars = c("pointer"), measure.vars = c("JN867756","JN867757",
                                                                                "MF621608","MG888045",
                                                                                "MH933760",
                                                                                "MN784123"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "JN867756"] <- "Homo sapiens cosavirus JN867756"
long.sim_nt$accession[long.sim_nt$accession == "JN867757"] <- "Homo sapiens cosavirus JN867757"
long.sim_nt$accession[long.sim_nt$accession == "MF621608"] <- "Homo sapiens cosavirus MF621608"
long.sim_nt$accession[long.sim_nt$accession == "MG888045"] <- "Miniopterus sp. mischivirus MG888045"
long.sim_nt$accession[long.sim_nt$accession == "MH933760"] <- "Miniopterus schreibersii shanbavirus MH933760"
long.sim_nt$accession[long.sim_nt$accession == "MN784123"] <- "Mandrillus leucophaeus cosavirus MN784123"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Homo sapiens cosavirus JN867756",
                                                                  "Homo sapiens cosavirus JN867757",
                                                                  "Homo sapiens cosavirus MF621608",
                                                                  "Miniopterus sp. mischivirus MG888045",
                                                                  "Miniopterus schreibersii shanbavirus MH933760",
                                                                  "Mandrillus leucophaeus cosavirus MN784123"))
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
mischi_bat_all_boot<-mischivirus_bat_all_boot/ictv_mischi+plot_layout(nrow=2,  heights = c(2, 0.30))
mischi_bat_all_boot

mischi_bat_all_boot<-as.ggplot(mischi_bat_all_boot)
mischi_bat_all_boot



#Sapelovirus
africa_sapelo_bootscan <- read.csv(file = "africa_sapelo_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_sapelo_bootscan)

#move to long
long.sim_nt <- melt(africa_sapelo_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","LC508226","NC_033820",
                                                                                       "OM104039"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_nt$accession[long.sim_nt$accession == "LC508226"] <- "Sus scrofa sapelovirus LC508226"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820"] <- "Eidolon helvum sapelovirus NC_033820"
long.sim_nt$accession[long.sim_nt$accession == "OM104039"] <- "Sus scrofa sapelovirus OM104039"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321",
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Sus scrofa sapelovirus LC508226",
                                                                  "Eidolon helvum sapelovirus NC_033820",
                                                                  "Sus scrofa sapelovirus OM104039"))
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
sapelo_bat_all_boot<-sapelovirus_bat_all_boot/ictv_sapelo_full+plot_layout(nrow=2,  heights = c(2, 0.30))
sapelo_bat_all_boot

sapelo_bat_all_boot<-as.ggplot(sapelo_bat_all_boot)
sapelo_bat_all_boot




#Sapovirus
africa_caliciviridae_bootscan <- read.csv(file = "africa_caliciviridae_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_caliciviridae_bootscan)

#move to long
long.sim_nt <- melt(africa_caliciviridae_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818340","OQ818345","OQ818347","OQ818348",
                                                                                    "JN699046", "KM056394","NC_033776","MW856468","MZ291974",
                                                                                    "OM105025"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818340"] <- "Eidolon helvum sapovirus OQ818340"
long.sim_nt$accession[long.sim_nt$accession == "OQ818345"] <- "Rousettus madagascariensis sapovirus OQ818345"
long.sim_nt$accession[long.sim_nt$accession == "OQ818347"] <- "Rousettus madagascariensis sapovirus OQ818347"
long.sim_nt$accession[long.sim_nt$accession == "OQ818348"] <- "Rousettus madagascariensis sapovirus OQ818348"
long.sim_nt$accession[long.sim_nt$accession == "JN699046"] <- "Homo sapiens norovirus JN699046"
long.sim_nt$accession[long.sim_nt$accession == "KM056394"] <- "Homo sapiens norovirus KM056394"
long.sim_nt$accession[long.sim_nt$accession == "NC_033776"] <- "Eidolon helvum sapovirus NC_033776"
long.sim_nt$accession[long.sim_nt$accession == "MW856468"] <- "Shellfish unclassified caliciviridae MW856468"
long.sim_nt$accession[long.sim_nt$accession == "MZ291974"] <- "Shellfish unclassified caliciviridae MZ291974"
long.sim_nt$accession[long.sim_nt$accession == "OM105025"] <- "Shellfish unclassified caliciviridae OM105025"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon helvum sapovirus OQ818340",
                                                                  "Rousettus madagascariensis sapovirus OQ818345",
                                                                  "Rousettus madagascariensis sapovirus OQ818347",
                                                                  "Rousettus madagascariensis sapovirus OQ818348",
                                                                  "Homo sapiens norovirus JN699046",
                                                                  "Homo sapiens norovirus KM056394", 
                                                                  "Eidolon helvum sapovirus NC_033776",
                                                                  "Shellfish unclassified caliciviridae MW856468",
                                                                  "Shellfish unclassified caliciviridae MZ291974",
                                                                  "Shellfish unclassified caliciviridae OM105025"))
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
sapo_bat_boot<-sapovirus_bat_boot/ictv_sapo_full+plot_layout(nrow=2,  heights = c(2, 0.30))
sapo_bat_boot

sapo_bat_boot<-as.ggplot(sapo_bat_boot)
sapo_bat_boot



#Teschovirus
africa_tescho_bootscan <- read.csv(file = "africa_tescho_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_tescho_bootscan)

#move to long
long.sim_nt <- melt(africa_tescho_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818318","OQ818324","OM966657","OM105029"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818318"] <- "Eidolon dupreanum teschovirus OQ818318"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324"] <- "Rousettus madagascariensis teschovirus OQ818324"
long.sim_nt$accession[long.sim_nt$accession == "OM966657"] <- "Sus scrofa teschovirus OM966657"
long.sim_nt$accession[long.sim_nt$accession == "OM105029"] <- "Sus scrofa teschovirus OM105029"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum teschovirus OQ818318",
                                                                  "Rousettus madagascariensis teschovirus OQ818324",
                                                                  "Sus scrofa teschovirus OM966657",
                                                                  "Sus scrofa teschovirus OM105029"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Rousettus madagascariensis teschovirus "), "OQ818323"))

teschovirus_bat_boot <- ggplot(long.sim_nt) + 
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

teschovirus_bat_boot

#put gene map with PySimPlot
tescho_boot<-teschovirus_bat_boot/ictv_tescho+plot_layout(nrow=2,  heights = c(2, 0.30))
tescho_boot

tescho_boot<-as.ggplot(tescho_boot)
tescho_boot







##Now put the whole figure together

#all bootscan together: all bat picornavirales
fig<-plot_grid(batpicorna_boot,
                           hep_bat_all_boot,
                           mischi_bat_all_boot,
                           sapelo_bat_all_boot,
                           sapo_bat_boot,
                           tescho_boot,
                    ncol=3,
                    labels="AUTO",  label_size = 23, align = "hv", axis="b")
fig

