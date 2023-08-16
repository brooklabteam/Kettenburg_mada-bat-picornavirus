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
bat <- read.csv("bat_genes.csv", header = T, stringsAsFactors = F)
bat$gene<-factor(bat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                        "VP1", "2A", "2B", "2C", "3A", "3B",
                                        "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                        "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                        "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the peptide files
bat_pep <- read.csv("bat_peptides.csv", header = T, stringsAsFactors = F)
bat_pep$gene<-factor(bat_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                "VP1", "2A", "2B", "2C", "3A", "3B",
                                                "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                                "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the feature file in case its needed
bat_feat <- read.csv("bat_features.csv", header = T, stringsAsFactors = F)
bat_feat$gene<-factor(bat_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                  "VP1", "2A", "2B", "2C", "3A", "3B",
                                                  "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                                  "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                  "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))

#Pick colors for genes
colz=c("5'UTR"="gold", "L"="royalblue","VP4"="paleturquoise3", "VP2"="skyblue1", "VP0"="royalblue4", "VP3"="steelblue1",
       "VP1"="cadetblue1", "2A"="palevioletred1", "2B"="red4", "2C"="palevioletred3", "3A"="tomato2", "3B"="plum",
       "3C"="rosybrown1", "3D"="pink2", 
       "Polyprotein"="azure3","Putative polyprotein"="mediumorchid1", "Non-structural polyprotein"="mediumorchid4", 
       "Putative minor structural protein" ="slateblue3", "Structural polyprotein"="lightgoldenrod1",
       "Hypothetical protein"="darkslategrey", "Similar to structural polyprotein"="mediumpurple1", "Similar to putative polyprotein"="mediumpurple3", 
       "Similar to polyprotein"="mediumpurple4","3'UTR"="yellow")


#Plot gene maps first

#Subset bat data by virus and dataset (all or full bat picornavirales)
bat_batpicornavirus<-subset(bat,molecule=="Bat picornavirus")
bat_hepatovirus<-subset(bat,molecule=="Hepatovirus")
bat_hepatovirus_full<-subset(bat_hepatovirus,type=="full")
bat_hepatovirus_all<-subset(bat_hepatovirus,type=="all")
bat_kunsagivirus<-subset(bat,molecule=="Kunsagivirus")
bat_mischivirus<-subset(bat,molecule=="Mischivirus")
bat_mischivirus_full<-subset(bat_mischivirus,type=="full")
bat_mischivirus_all<-subset(bat_mischivirus,type=="all")
bat_sapelovirus<-subset(bat,molecule=="Sapelovirus")
bat_sapelovirus_full<-subset(bat_sapelovirus,type=="full")
bat_sapelovirus_all<-subset(bat_sapelovirus,type=="all")
bat_teschovirus<-subset(bat,molecule=="Teschovirus")
bat_sapovirus<-subset(bat,molecule=="Sapovirus")
bat_felispicornalike<-subset(bat,molecule=="Felisavirus and Picorna-like virus")
bat_felisavirus<-subset(bat,molecule=="Felisavirus")
bat_picornalike<-subset(bat,molecule=="Picorna-like virus")


bat_batpicornavirus_pep<-subset(bat_pep,molecule=="Bat picornavirus")
bat_hepatovirus_pep<-subset(bat_pep,molecule=="Hepatovirus")
bat_hepatovirus_full_pep<-subset(bat_hepatovirus_pep,type=="full")
bat_hepatovirus_all_pep<-subset(bat_hepatovirus_pep,type=="all")
bat_kunsagivirus_pep<-subset(bat_pep,molecule=="Kunsagivirus")
bat_mischivirus_pep<-subset(bat_pep,molecule=="Mischivirus")
bat_mischivirus_full_pep<-subset(bat_mischivirus_pep,type=="full")
bat_mischivirus_all_pep<-subset(bat_mischivirus_pep,type=="all")
bat_sapelovirus_pep<-subset(bat_pep,molecule=="Sapelovirus")
bat_sapelovirus_full_pep<-subset(bat_sapelovirus_pep,type=="full")
bat_sapelovirus_all_pep<-subset(bat_sapelovirus_pep,type=="all")
bat_teschovirus_pep<-subset(bat_pep,molecule=="Teschovirus")


bat_batpicornavirus_feat<-subset(bat_feat,molecule=="Bat picornavirus")
bat_hepatovirus_feat<-subset(bat_feat,molecule=="Hepatovirus")
bat_hepatovirus_full_feat<-subset(bat_hepatovirus_feat,type=="full")
bat_hepatovirus_all_feat<-subset(bat_hepatovirus_feat,type=="all")
bat_kunsagivirus_feat<-subset(bat_feat,molecule=="Kunsagivirus")
bat_mischivirus_feat<-subset(bat_feat,molecule=="Mischivirus")
bat_mischivirus_full_feat<-subset(bat_mischivirus_feat,type=="full")
bat_mischivirus_all_feat<-subset(bat_mischivirus_feat,type=="all")
bat_sapelovirus_feat<-subset(bat_feat,molecule=="Sapelovirus")
bat_sapelovirus_full_feat<-subset(bat_sapelovirus_feat,type=="full")
bat_sapelovirus_all_feat<-subset(bat_sapelovirus_feat,type=="all")
bat_teschovirus_feat<-subset(bat_feat,molecule=="Teschovirus")
bat_sapovirus_feat<-subset(bat_feat,molecule=="Sapovirus")
bat_felispicornalike_feat<-subset(bat_feat,molecule=="Felisavirus and Picorna-like virus")
bat_felisavirus_feat<-subset(bat_feat,molecule=="Felisavirus")
bat_picornalike_feat<-subset(bat_feat,molecule=="Picorna-like virus")



##Now the bat gene maps



#all seq over 3kb
bat_batpicorna_all<-ggplot(bat_batpicornavirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_batpicornavirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_batpicornavirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1300,7940),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_batpicorna_all


felispicorna_all<-ggplot(bat_felispicornalike, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_text(data=bat_felispicornalike_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(3600,11970),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
felispicorna_all


bat_hepato_all<-ggplot(bat_hepatovirus_all, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_hepatovirus_all_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                  xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_hepatovirus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(800,6362),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_hepato_all


bat_mischi_all<-ggplot(bat_mischivirus_all, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_mischivirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_mischivirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_mischivirus_all_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                             xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_mischivirus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1590,8850),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_mischi_all


bat_sapelo_all<-ggplot(bat_sapelovirus_all, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_sapelovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_sapelovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_sapelovirus_all_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                  xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_sapelovirus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1290,8200),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_sapelo_all


bat_sapo_all<-ggplot(bat_sapovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_sapovirus_full_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_sapovirus_full_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_text(data=bat_sapovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1050,7530),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_sapo_all


bat_tescho_all<-ggplot(bat_teschovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_teschovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_teschovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_teschovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                             xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_teschovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(4350,7400),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_tescho_all



##now do all full reference genome maps
bat_felis_full<-ggplot(bat_felisavirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_text(data=bat_felisavirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(540,8970),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_felis_full


bat_hepato_full<-ggplot(bat_hepatovirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_hepatovirus_full_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_hepatovirus_full_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(160,6175),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_hepato_full


bat_kun_full<-ggplot(bat_kunsagivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_mischivirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_mischivirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_kunsagivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_kunsagivirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(490,8780),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_kun_full


bat_mischi_full<-ggplot(bat_mischivirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_sapelovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_sapelovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_mischivirus_full_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_mischivirus_full_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(410,8870),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_mischi_full


bat_picornalike_full<-ggplot(bat_picornalike, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_sapovirus_full_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_sapovirus_full_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_text(data=bat_picornalike_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(530,8860),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_picornalike_full


bat_sapelo_full<-ggplot(bat_sapelovirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_teschovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_teschovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_sapelovirus_full_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                             xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_sapelovirus_full_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(110,8060),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_sapelo_full





##Now get the plots from bootscan
setwd("~/Desktop/bootscan/output")


#Bat picornavirus 
bat_picorna_bat_boot <- read.csv(file = "bat_picorna_bat_bootscan.csv", header = T, stringsAsFactors = F)
head(bat_picorna_bat_boot)

#move to long
long.sim_nt <- melt(bat_picorna_bat_boot, id.vars = c("pointer"), measure.vars = c("OQ818328","HQ595340","HQ595342","HQ595344"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818328"] <- "Rousettus madagascariensis picornavirus OQ818328"
long.sim_nt$accession[long.sim_nt$accession == "HQ595340"] <- "Bat picornavirus 1 HQ595340"
long.sim_nt$accession[long.sim_nt$accession == "HQ595342"] <- "Bat picornavirus 2 HQ595342"
long.sim_nt$accession[long.sim_nt$accession == "HQ595344"] <- "Bat picornavirus 3 HQ595344"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rousettus madagascariensis picornavirus OQ818328", "Bat picornavirus 1 HQ595340",
                                                                  "Bat picornavirus 2 HQ595342","Bat picornavirus 3 HQ595344"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("Reference: ",italic("Rousettus madagascariensis picornavirus "), "OQ818325"))

batpicorna_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "top", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold"), ) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
  #                    labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_bat_boot

#put gene map with PySimPlot
batpicorna_boot<-batpicorna_bat_boot/bat_batpicorna+plot_layout(nrow=2,  heights = c(2, 0.30))
batpicorna_boot

batpicorna_boot<-as.ggplot(batpicorna_boot)
batpicorna_boot




#Cheravirus RNA2
cheravirus_bat_boot <- read.csv(file = "chera_bat_bootscan.csv", header = T, stringsAsFactors = F) 
head(cheravirus_bat_boot)

#move to long
long.sim_nt <- melt(cheravirus_bat_boot, id.vars = c("pointer"), measure.vars = c("AB030941","AJ621358","KT692953","DQ143875","MK153132","MW582786"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "AB030941"] <- "Cheravirus mali AB030941"
long.sim_nt$accession[long.sim_nt$accession == "AJ621358"] <- "Cheravirus avii AJ621358"
long.sim_nt$accession[long.sim_nt$accession == "KT692953"] <- "Cheravirus ribis KT692953"
long.sim_nt$accession[long.sim_nt$accession == "DQ143875"] <- "Cheravirus pruni DQ143875"
long.sim_nt$accession[long.sim_nt$accession == "MK153132"] <- "Arracacha virus B MK153132"
long.sim_nt$accession[long.sim_nt$accession == "MW582786"] <- "Cheravirus arracaciae MW582786"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Cheravirus mali AB030941", "Cheravirus avii AJ621358",
                                                                  "Cheravirus ribis KT692953","Cheravirus pruni DQ143875",
                                                                  "Arracacha virus B MK153132","Cheravirus arracaciae MW582786"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100



## Nucleotide
title<-expression(paste("Reference: ",italic("Pteropus rufus cheravirus "), "RNA2 OQ818330"))

cheravirus_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "top", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,1000/3.055,2000/3.055,3000/3.055), 
  #                    labels = c(0,1000, 2000,3000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

cheravirus_bat_boot

#put gene map with PySimPlot
chera_boot<-cheravirus_bat_boot/bat_chera+plot_layout(nrow=2,  heights = c(2, 0.30))
chera_boot

chera_boot<-as.ggplot(chera_boot)
chera_boot






#Hepatovirus
hepato_bat_boot <- read.csv(file = "hepato_bat_bootscan.csv", header = T, stringsAsFactors = F) #animo acid
head(hepato_bat_boot)

#move to long
long.sim_nt <- melt(hepato_bat_boot, id.vars = c("pointer"), measure.vars = c("HPA","KR703607","KT452637","KT452658",
                                                                                     "KT452685","KT452714","KT452730", "KT452735",
                                                                                     "KT452742"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "HPA"] <- "Hepatovirus A HPA"
long.sim_nt$accession[long.sim_nt$accession == "KR703607"] <- "Hepatovirus B KR703607"
long.sim_nt$accession[long.sim_nt$accession == "KT452637"] <- "Hepatovirus D KT452637"
long.sim_nt$accession[long.sim_nt$accession == "KT452658"] <- "Hepatovirus I KT452658"
long.sim_nt$accession[long.sim_nt$accession == "KT452685"] <- "Hepatovirus F KT452685"
long.sim_nt$accession[long.sim_nt$accession == "KT452714"] <- "Hepatovirus H KT452714"
long.sim_nt$accession[long.sim_nt$accession == "KT452730"] <- "Hepatovirus G KT452730"
long.sim_nt$accession[long.sim_nt$accession == "KT452735"] <- "Hepatovirus E KT452735"
long.sim_nt$accession[long.sim_nt$accession == "KT452742"] <- "Hepatovirus C KT452742"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Hepatovirus A HPA", "Hepatovirus B KR703607",
                                                                  "Hepatovirus C KT452742","Hepatovirus D KT452637",
                                                                  "Hepatovirus E KT452735","Hepatovirus F KT452685",
                                                                  "Hepatovirus G KT452730", "Hepatovirus H KT452714",
                                                                  "Hepatovirus I KT452658"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum hepatovirus "), "OQ818337"))

hepatovirus_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
  annotate("point",x=1220, y=0.30, size=3, color="black")+
  annotate("point",x=1300, y=0.45, size=3, color="black")+
  annotate("point",x=2090, y=0.35, size=3, color="black")+
  annotate("point",x=5270, y=0.30, size=3, color="black")+
  # geom_vline(xintercept=1220, linetype="dashed")+
  # geom_vline(xintercept=1300, linetype="dashed")+
  # geom_vline(xintercept=2090, linetype="dashed")+
  # geom_vline(xintercept=5270, linetype="dashed")+
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "top", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_bat_boot


#put gene map with PySimPlot
hep_bat_boot<-hepatovirus_bat_boot/bat_hepato+plot_layout(nrow=2,  heights = c(2, 0.30))
hep_bat_boot

hep_bat_boot<-as.ggplot(hep_bat_boot)
hep_bat_boot



#kobuvirus
kobuvirus_bat_boot <- read.csv(file = "kobu_bat_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(kobuvirus_bat_boot)

#move to long
long.sim_nt <- melt(kobuvirus_bat_boot, id.vars = c("pointer"), measure.vars = c("AB040749","AB084788",
                                                                                "EU787450","KJ641686",
                                                                                "LC055961",
                                                                                "OP287812"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "AB040749"] <- "Aichivirus A AB040749"
long.sim_nt$accession[long.sim_nt$accession == "AB084788"] <- "Aichivirus B AB084788"
long.sim_nt$accession[long.sim_nt$accession == "EU787450"] <- "Aichivirus C EU787450"
long.sim_nt$accession[long.sim_nt$accession == "KJ641686"] <- "Aichivirus F KJ641686"
long.sim_nt$accession[long.sim_nt$accession == "LC055961"] <- "Aichivirus D LC055961"
long.sim_nt$accession[long.sim_nt$accession == "OP287812"] <- "Eidolon dupreanum kobuvirus OP287812"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Aichivirus A AB040749", "Aichivirus B AB084788",
                                                                  "Aichivirus C EU787450","Aichivirus D LC055961",
                                                                  "Aichivirus F KJ641686",
                                                                  "Eidolon dupreanum kobuvirus OP287812"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum kobuvirus "), "OQ818322"))

kobuvirus_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "top", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobuvirus_bat_boot

#put gene map with PySimPlot
kobu_bat_boot<-kobuvirus_bat_boot/bat_kobu+plot_layout(nrow=2,  heights = c(2, 0.30))
kobu_bat_boot

kobu_bat_bppt<-as.ggplot(kobu_bat_boot)
kobu_bat_boot





#kunsagivirus
kunsagivirus_bat_boot <- read.csv(file = "kun_bat_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(kunsagivirus_bat_boot)

#move to long
long.sim_nt <- melt(kunsagivirus_bat_boot, id.vars = c("pointer"), measure.vars = c("KC935379","KX644936",
                                                                                   "KY670597"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "KC935379"] <- "Kunsagivirus A KC935379"
long.sim_nt$accession[long.sim_nt$accession == "KX644936"] <- "Kunsagivirus B KX644936"
long.sim_nt$accession[long.sim_nt$accession == "KY670597"] <- "Kunsagivirus C KY670597"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Kunsagivirus A KC935379", "Kunsagivirus B KX644936",
                                                                  "Kunsagivirus C KY670597"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum kunsagivirus "), "OQ818317"))

kunsagivirus_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
  annotate("point",x=3230, y=0.57, size=3, color="black")+
  #geom_vline(xintercept=3230, linetype="dashed")+
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("%permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "top", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagivirus_bat_boot


#put gene map with PySimPlot
kun_bat_boot<-kunsagivirus_bat_boot/bat_kun+plot_layout(nrow=2,  heights = c(2, 0.30))
kun_bat_boot

kun_bat_boot<-as.ggplot(kun_bat_boot)
kun_bat_boot




#Mischivirus
mischivirus_bat_boot <- read.csv(file = "mischi_bat_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(mischivirus_bat_boot)

#move to long
long.sim_nt <- melt(mischivirus_bat_boot, id.vars = c("pointer"), measure.vars = c("JQ814851","KP054273",
                                                                                  "KP100644","KY512802",
                                                                                  "MF352410"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "JQ814851"] <- "Mischivirus A JQ814851"
long.sim_nt$accession[long.sim_nt$accession == "KP054273"] <- "Mischivirus B KP054273"
long.sim_nt$accession[long.sim_nt$accession == "KP100644"] <- "Mischivirus C KP100644"
long.sim_nt$accession[long.sim_nt$accession == "KY512802"] <- "Mischivirus D KY512802"
long.sim_nt$accession[long.sim_nt$accession == "MF352410"] <- "Mischivirus E MF352410"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Mischivirus A JQ814851", "Mischivirus B KP054273",
                                                                  "Mischivirus C KP100644", "Mischivirus D KY512802",
                                                                  "Mischivirus E MF352410"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Pteropus rufus mischivirus "), "OQ818316"))

mischivirus_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "top", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_bat_boot

#put gene map with PySimPlot
mischi_bat_boot<-mischivirus_bat_boot/bat_mischi+plot_layout(nrow=2,  heights = c(2, 0.30))
mischi_bat_boot

mischi_bat_boot<-as.ggplot(mischi_bat_boot)
mischi_bat_boot



#Sapelovirus full
sapelovirus_bat_boot <- read.csv(file = "sapelo_bat_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(sapelovirus_bat_boot)

#move to long
long.sim_nt <- melt(sapelovirus_bat_boot, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","AF406813","AY064708",
                                                                                       "NC_033820"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_nt$accession[long.sim_nt$accession == "AF406813"] <- "Sapelovirus A AF406813"
long.sim_nt$accession[long.sim_nt$accession == "AY064708"] <- "Sapelovirus B AY064708"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820"] <- "Eidolon helvum sapelovirus NC_033820"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321",
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Sapelovirus A AF406813","Sapelovirus B AY064708",
                                                                  "Eidolon helvum sapelovirus NC_033820"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapelovirus "), "OQ818320"))

sapelovirus_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "top", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_bat_boot

#put gene map with PySimPlot
sapelo_bat_boot<-sapelovirus_bat_boot/bat_sapelo_full+plot_layout(nrow=2,  heights = c(2, 0.30))
sapelo_bat_boot

sapelo_bat_boot<-as.ggplot(sapelo_bat_boot)
sapelo_bat_boot





#Sapovirus
sapovirus_boot <- read.csv(file = "sapo_all_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(sapovirus_boot)

#move to long
long.sim_nt <- melt(sapovirus_boot, id.vars = c("pointer"), measure.vars = c("HM002617","KX759618.1","KX759621.1","KX759622.1",
                                                                                    "NC_033776.1", "KX759623.1"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "HM002617"] <- "Sapporo virus HM002617"
long.sim_nt$accession[long.sim_nt$accession == "KX759618.1"] <- "Eidolon helvum sapovirus KX759618"
long.sim_nt$accession[long.sim_nt$accession == "KX759621.1"] <- "Eidolon helvum sapovirus KX759621"
long.sim_nt$accession[long.sim_nt$accession == "KX759622.1"] <- "Eidolon helvum sapovirus KX759622"
long.sim_nt$accession[long.sim_nt$accession == "NC_033776.1"] <- "Eidolon helvum sapovirus NC_033776"
long.sim_nt$accession[long.sim_nt$accession == "KX759623.1"] <- "Eidolon helvum sapovirus KX759623"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Sapporo virus HM002617",
                                                                  "Eidolon helvum sapovirus KX759618",
                                                                  "Eidolon helvum sapovirus KX759621",
                                                                  "Eidolon helvum sapovirus KX759622",
                                                                  "Eidolon helvum sapovirus KX759623", 
                                                                  "Eidolon helvum sapovirus NC_033776"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapovirus "), "OQ818319"))

sapovirus_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  annotate("point",x=2410, y=0.45, size=3, color="black")+
  annotate("point",x=2520, y=0.42, size=3, color="black")+
  annotate("point",x=2820, y=0.46, size=3, color="black")+
  annotate("point",x=3740, y=0.47, size=3, color="black")+
  annotate("point",x=4600, y=0.49, size=3, color="black")+
  annotate("point",x=7180, y=0.4, size=3, color="black")+
  annotate("point",x=7240, y=0.38, size=3, color="black")+
  geom_hline(yintercept=0.30, linetype="dashed", color="lightgrey")+
  #geom_vline(xintercept=2410, linetype="dashed")+
  #geom_vline(xintercept=2520, linetype="dashed")+
  #geom_vline(xintercept=2820, linetype="dashed")+
  #geom_vline(xintercept=3740, linetype="dashed")+
  #geom_vline(xintercept=4600, linetype="dashed")+
  #geom_vline(xintercept=7180, linetype="dashed")+
  #geom_vline(xintercept=7240, linetype="dashed")+
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "top", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_bat_boot

#put gene map with PySimPlot
sapo_bat_boot<-sapovirus_bat_boot/bat_sapo_full+plot_layout(nrow=2,  heights = c(2, 0.30))
sapo_bat_boot

sapo_bat_boot<-as.ggplot(sapo_bat_boot)
sapo_bat_boot





#Teschovirus
teschovirus_bat_boot <- read.csv(file = "tescho_bat_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(teschovirus_bat_boot)

#move to long
long.sim_nt <- melt(teschovirus_bat_boot, id.vars = c("pointer"), measure.vars = c("OQ818323","OQ818324","LC386158","MG875515","MT295502"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818323"] <- "Rousettus madagascariensis teschovirus OQ818323"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324"] <- "Rousettus madagascariensis teschovirus OQ818324"
long.sim_nt$accession[long.sim_nt$accession == "LC386158"] <- "Teschovirus A LC386158"
long.sim_nt$accession[long.sim_nt$accession == "MG875515"] <- "Teschovirus B MG875515"
long.sim_nt$accession[long.sim_nt$accession == "MT295502"] <- "Teschovirus A6 MT295502"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rousettus madagascariensis teschovirus OQ818323",
                                                                  "Rousettus madagascariensis teschovirus OQ818324",
                                                                  "Teschovirus A LC386158",
                                                                  "Teschovirus B MG875515",
                                                                  "Teschovirus A6 MT295502"))
#and plot
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
        legend.position = "top", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 7, ),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

teschovirus_bat_boot

#put gene map with PySimPlot
tescho_boot<-teschovirus_bat_boot/bat_tescho+plot_layout(nrow=2,  heights = c(2, 0.30))
tescho_boot

tescho_boot<-as.ggplot(tescho_boot)
tescho_boot




##Now put the whole figure together

#all bootscan together
bootscan<-plot_grid(batpicorna_boot,
                    mischi_bat_boot,
                    chera_boot,
                    hep_bat_boot,
                    kobu_bat_boot,
                    kun_bat_boot,
                    sapelo_bat_boot, 
                    sapo_bat_boot,
                    tescho_boot,
                    ncol=3,
                    labels="AUTO")
bootscan


#in Fig 3
bootscan_fig<-plot_grid(hep_bat_boot,
                        kun_bat_boot, 
                        sapo_bat_boot,
                        nrow=1,
                        labels="AUTO")
bootscan_fig


#excluded from Fig 3 
bootscan_supp<-plot_grid(batpicorna_boot,
                         chera_boot,
                         mischi_bat_boot,
                         kobu_bat_boot,
                         sapelo_bat_boot,
                         tescho_boot,
                         ncol=3,
                         labels="AUTO")
bootscan_supp

