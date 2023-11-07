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
ictv <- read.csv("ictv_blast_genes_trimmed.csv", header = T, stringsAsFactors = F)
ictv$gene<-factor(ictv$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                        "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                        "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Minor structural protein", "Non-structural polyprotein",
                                        "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                        "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the peptide files
ictv_pep <- read.csv("ictv_blast_peptides_trimmed.csv", header = T, stringsAsFactors = F)
ictv_pep$gene<-factor(ictv_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                                "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Minor structural protein", "Non-structural polyprotein",
                                                "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the feature file in case its needed
ictv_feat <- read.csv("ictv_blast_features_trimmed.csv", header = T, stringsAsFactors = F)
ictv_feat$gene<-factor(ictv_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                  "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                                  "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Minor structural protein", "Non-structural polyprotein",
                                                  "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                  "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))

#Pick colors for genes
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
ictv_batpicornavirus_full<-subset(ictv_batpicornavirus,type=="full")
ictv_batpicornavirus_all<-subset(ictv_batpicornavirus,type=="all")
ictv_hepatovirus<-subset(ictv,molecule=="Hepatovirus")
ictv_kobuvirus<-subset(ictv,molecule=="Kobuvirus")
ictv_kunsagivirus<-subset(ictv,molecule=="Kunsagivirus")
ictv_mischivirus<-subset(ictv,molecule=="Mischivirus")
ictv_sapovirus<-subset(ictv,molecule=="Sapovirus")
ictv_sapovirus_full<-subset(ictv_sapovirus,type=="full")
ictv_sapovirus_all<-subset(ictv_sapovirus,type=="all")
ictv_sapelovirus<-subset(ictv,molecule=="Sapelovirus")
ictv_sapelovirus_full<-subset(ictv_sapelovirus,type=="full")
ictv_sapelovirus_p1<-subset(ictv_sapelovirus,type=="p1")
ictv_sapelovirus_p2<-subset(ictv_sapelovirus,type=="p2")
ictv_teschovirus<-subset(ictv,molecule=="Teschovirus")

ictv_batpicornavirus_pep<-subset(ictv_pep,molecule=="Bat picornavirus")
ictv_batpicornavirus_full_pep<-subset(ictv_batpicornavirus_pep,type=="full")
ictv_batpicornavirus_all_pep<-subset(ictv_batpicornavirus_pep,type=="all")
ictv_hepatovirus_pep<-subset(ictv_pep, molecule=="Hepatovirus")
ictv_kobuvirus_pep<-subset(ictv_pep,molecule=="Kobuvirus")
ictv_kunsagivirus_pep<-subset(ictv_pep,molecule=="Kunsagivirus")
ictv_mischivirus_pep<-subset(ictv_pep,molecule=="Mischivirus")
ictv_sapovirus_pep<-subset(ictv_pep, molecule=="Sapovirus")
ictv_sapovirus_full_pep<-subset(ictv_sapovirus_pep,type=="full")
ictv_sapovirus_all_pep<-subset(ictv_sapovirus_pep,type=="all")
ictv_sapelovirus_pep<-subset(ictv_pep,molecule=="Sapelovirus")
ictv_sapelovirus_pep_full<-subset(ictv_sapelovirus_pep,type=="full")
ictv_sapelovirus_pep_p1<-subset(ictv_sapelovirus_pep,type=="p1")
ictv_sapelovirus_pep_p2<-subset(ictv_sapelovirus_pep,type=="p2")
ictv_teschovirus_pep<-subset(ictv_pep,molecule=="Teschovirus")

ictv_batpicornavirus_feat<-subset(ictv_feat,molecule=="Bat picornavirus")
ictv_batpicornavirus_full_feat<-subset(ictv_batpicornavirus_feat,type=="full")
ictv_batpicornavirus_all_feat<-subset(ictv_batpicornavirus_feat,type=="all")
ictv_hepatovirus_feat<-subset(ictv_feat,molecule=="Hepatovirus")
ictv_kobuvirus_feat<-subset(ictv_feat,molecule=="Kobuvirus")
ictv_kunsagivirus_feat<-subset(ictv_feat,molecule=="Kunsagivirus")
ictv_mischivirus_feat<-subset(ictv_feat,molecule=="Mischivirus")
ictv_sapovirus_feat<-subset(ictv_feat,molecule=="Sapovirus")
ictv_sapovirus_full_feat<-subset(ictv_sapovirus_feat,type=="full")
ictv_sapovirus_all_feat<-subset(ictv_sapovirus_feat,type=="all")
ictv_sapelovirus_feat<-subset(ictv_feat,molecule=="Sapelovirus")
ictv_sapelovirus_full_feat<-subset(ictv_sapelovirus_feat,type=="full")
ictv_sapelovirus_feat_p1<-subset(ictv_sapelovirus_feat,type=="p1")
ictv_sapelovirus_feat_p2<-subset(ictv_sapelovirus_feat,type=="p2")
ictv_teschovirus_feat<-subset(ictv_feat,molecule=="Teschovirus")

#plot ictv and blast plots
ictv_batpicorna_all<-ggplot(ictv_batpicornavirus_all, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_batpicornavirus_all_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                  xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_batpicornavirus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(4545,5605),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_batpicorna_all


ictv_batpicorna_full<-ggplot(ictv_batpicornavirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_batpicornavirus_full_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                      xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_batpicornavirus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1,7700),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_batpicorna_full


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
  scale_x_continuous(limits=c(0,6250),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,8650),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,7900),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapelo_full


ictv_sapelo_p1<-ggplot(ictv_sapelovirus_p1, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_sapelovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_sapelovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_sapelovirus_pep_p1, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                   xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_sapelovirus_feat_p1,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(810,1975),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapelo_p1


ictv_sapelo_p2<-ggplot(ictv_sapelovirus_p2, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_sapelovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_sapelovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_sapelovirus_pep_p2, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_sapelovirus_feat_p2,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(5030,5860),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapelo_p2



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
  scale_x_continuous(limits=c(0,7450),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_tescho


ictv_sapo_all<-ggplot(ictv_sapovirus_all, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_sapovirus_full_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_sapovirus_full_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_sapovirus_all_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_sapovirus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(2700,4290),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapo_all


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
  scale_x_continuous(limits=c(700,7950),expand=c(0,0))+
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

colzpalette<-c("coral1","cyan3","hotpink1","gold","orange3","dodgerblue2","firebrick1","darkorchid3","skyblue1")

#Bat picornavirus all
africa_batpicorna_all_bootscan <- read.csv(file = "africa_batpicorna_all_bootscan.csv", header = T, stringsAsFactors = F)
head(africa_batpicorna_all_bootscan)

#move to long
long.sim_nt <- melt(africa_batpicorna_all_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818325","OQ818346","JX437642","KF040078"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818325"] <- "R. madagascariensis roupivirus MIZ214"
long.sim_nt$accession[long.sim_nt$accession == "OQ818346"] <- "R. madagascariensis roupivirus MIZ194"
long.sim_nt$accession[long.sim_nt$accession == "JX437642"] <- "Homo sapiens enterovirus: JX437642"
long.sim_nt$accession[long.sim_nt$accession == "KF040078"] <- "Pan troglodytes enterovirus: KF040078"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. madagascariensis roupivirus MIZ214", 
                                                                  "R. madagascariensis roupivirus MIZ194",
                                                                  "Homo sapiens enterovirus: JX437642","Pan troglodytes enterovirus: KF040078"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("Reference: . madagascariensis roupivirus MIZ240"))

batpicorna_africa_all_boot <- ggplot(long.sim_nt) + 
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

batpicorna_africa_all_boot

#put gene map with PySimPlot
batpicorna_all_boot<-batpicorna_africa_all_boot/ictv_batpicorna_all+plot_layout(nrow=2,  heights = c(1, 0.30))
batpicorna_all_boot

batpicorna_all_boot<-as.ggplot(batpicorna_all_boot)
batpicorna_all_boot



#Bat picornavirus full
africa_batpicorna_full_bootscan <- read.csv(file = "africa_batpicorna_full_bootscan.csv", header = T, stringsAsFactors = F)
head(africa_batpicorna_full_bootscan)

#move to long
long.sim_nt <- melt(africa_batpicorna_full_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818325","JX437642","KF040078"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818325"] <- "R. madagascariensis roupivirus MIZ214"
long.sim_nt$accession[long.sim_nt$accession == "JX437642"] <- "Homo sapiens enterovirus: JX437642"
long.sim_nt$accession[long.sim_nt$accession == "KF040078"] <- "Pan troglodytes enterovirus: KF040078"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. madagascariensis roupivirus MIZ214", 
                                                                  "Homo sapiens enterovirus: JX437642",
                                                                  "Pan troglodytes enterovirus: KF040078"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("R. madagascariensis roupivirus MIZ240"))

batpicorna_africa_full_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
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
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
  #                    labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_africa_full_boot

#put gene map with PySimPlot
batpicorna_full_boot<-batpicorna_africa_full_boot/ictv_batpicorna_full+plot_layout(nrow=2,  heights = c(1, 0.30))
batpicorna_full_boot

batpicorna_full_boot<-as.ggplot(batpicorna_full_boot)
batpicorna_full_boot



#Hepatovirus
africa_hepato_bootscan <- read.csv(file = "africa_hepato_bootscan.csv", header = T, stringsAsFactors = F) #animo acid
head(africa_hepato_bootscan)

#move to long
long.sim_nt <- melt(africa_hepato_bootscan, id.vars = c("pointer"), measure.vars = c("KT452729",
                                                                                     "NC_028366","NC_038313","NC_038314"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "KT452729"] <- "Rhinolophus landeri hepatovirus: KT452729"
long.sim_nt$accession[long.sim_nt$accession == "NC_028366"] <- "Eidolon helvum hepatovirus: NC_028366"
long.sim_nt$accession[long.sim_nt$accession == "NC_038313"] <- "Miniopterus sp. hepatovirus: NC_038313"
long.sim_nt$accession[long.sim_nt$accession == "NC_038314"] <- "Lophuromys sikapusi hepatovirus: NC_038314"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rhinolophus landeri hepatovirus: KT452729",
                                                                  "Eidolon helvum hepatovirus: NC_028366",
                                                                  "Miniopterus sp. hepatovirus: NC_038313",
                                                                  "Lophuromys sikapusi hepatovirus: NC_038314"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: E. dupreanum hepatovirus KEL148"))

hepatovirus_bat_all_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
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
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_bat_all_boot


#put gene map with PySimPlot
hep_bat_all_boot<-hepatovirus_bat_all_boot/ictv_hepato+plot_layout(nrow=2,  heights = c(1, 0.30))
hep_bat_all_boot

hep_bat_all_boot<-as.ggplot(hep_bat_all_boot)
hep_bat_all_boot




#kunsagivirus
africa_kun_bootscan <- read.csv(file = "africa_kun_bootscan.csv", header = T, stringsAsFactors = F) #animo acid
head(africa_kun_bootscan)

#move to long
long.sim_nt <- melt(africa_kun_bootscan, id.vars = c("pointer"), measure.vars = c("NC_033818",
                                                                                     "NC_034206","HF677705"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_033818"] <- "Eidolon helvum kunsagivirus: NC_033818"
long.sim_nt$accession[long.sim_nt$accession == "NC_034206"] <- "Papio cynocephalus kunsagivirus: NC_034206"
long.sim_nt$accession[long.sim_nt$accession == "HF677705"] <- "Hylomyscus parechovirus: HF677705"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon helvum kunsagivirus: NC_033818",
                                                                  "Papio cynocephalus kunsagivirus: NC_034206",
                                                                  "Hylomyscus parechovirus: HF677705"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: E. dupreanum kunsagivirus KEL148"))

kunsagivirus_bat_all_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
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
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagivirus_bat_all_boot


#put gene map with PySimPlot
kun_bat_all_boot<-kunsagivirus_bat_all_boot/ictv_kun+plot_layout(nrow=2,  heights = c(1, 0.30))
kun_bat_all_boot

kun_bat_all_boot<-as.ggplot(kun_bat_all_boot)
kun_bat_all_boot



#mischivirus
africa_mischi_bootscan <- read.csv(file = "africa_mischi_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_mischi_bootscan)

#move to long
long.sim_nt <- melt(africa_mischi_bootscan, id.vars = c("pointer"), measure.vars = c("JN867757",
                                                                                "MG888045",
                                                                                "MN784123"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "JN867757"] <- "Homo sapiens cosavirus: JN867757"
long.sim_nt$accession[long.sim_nt$accession == "MG888045"] <- "Miniopterus sp. mischivirus: MG888045"
long.sim_nt$accession[long.sim_nt$accession == "MN784123"] <- "Mandrillus leucophaeus cosavirus: MN784123"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Homo sapiens cosavirus: JN867757",
                                                                  "Miniopterus sp. mischivirus: MG888045",
                                                                  "Mandrillus leucophaeus cosavirus: MN784123"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: P. rufus mischivirus AMB150"))

mischivirus_bat_all_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
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
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_bat_all_boot

#put gene map with PySimPlot
mischi_bat_all_boot<-mischivirus_bat_all_boot/ictv_mischi+plot_layout(nrow=2,  heights = c(1, 0.30))
mischi_bat_all_boot

mischi_bat_all_boot<-as.ggplot(mischi_bat_all_boot)

mischi_bat_all_boot



#Sapelovirus
africa_sapelo_full_bootscan <- read.csv(file = "africa_sapelo_full_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_sapelo_full_bootscan)

#move to long
long.sim_nt <- melt(africa_sapelo_full_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","NC_033820"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "E. dupreanum sapelovirus KEL272"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "R. madagascariensis sapelovirus MIZ243"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820"] <- "Eidolon helvum sapelovirus: NC_033820"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum sapelovirus KEL272",
                                                                  "R. madagascariensis sapelovirus MIZ243",
                                                                  "Eidolon helvum sapelovirus: NC_033820"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: E. dupreanum sapelovirus KEL233"))

sapelovirus_bat_full_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
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
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_bat_full_boot

#put gene map with PySimPlot
sapelo_bat_full_boot<-sapelovirus_bat_full_boot/ictv_sapelo_full+plot_layout(nrow=2,  heights = c(1, 0.30))
sapelo_bat_full_boot

sapelo_bat_full_boot<-as.ggplot(sapelo_bat_full_boot)
sapelo_bat_full_boot


#Sapelovirus p2
africa_sapelo_p1_bootscan <- read.csv(file = "africa_sapelo_p1_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_sapelo_p1_bootscan)

#move to long
long.sim_nt <- melt(africa_sapelo_p1_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","NC_033820",
                                                                                        "OQ818343","LC508226","OM104039"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "E. dupreanum sapelovirus KEL272"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "R. madagascariensis sapelovirus MIZ243"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820"] <- "Eidolon helvum sapelovirus: NC_033820"
long.sim_nt$accession[long.sim_nt$accession == "OQ818343"] <- "E. dupreanum sapelovirus KEL275B"
long.sim_nt$accession[long.sim_nt$accession == "LC508226"] <- "Porcine sapelovirus: LC508226"
long.sim_nt$accession[long.sim_nt$accession == "OM104039"] <- "Porcine sapelovirus: OM104039"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum sapelovirus KEL272",
                                                                  "R. madagascariensis sapelovirus MIZ243",
                                                                  "Eidolon helvum sapelovirus: NC_033820",
                                                                  "E. dupreanum sapelovirus KEL275B",
                                                                  "Porcine sapelovirus: LC508226",
                                                                  "Porcine sapelovirus: OM104039"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: E. dupreanum sapelovirus KEL233"))

sapelovirus_bat_p1_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
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

sapelovirus_bat_p1_boot

#put gene map with PySimPlot
sapelo_bat_p1_boot<-sapelovirus_bat_p1_boot/ictv_sapelo_p1+plot_layout(nrow=2,  heights = c(1, 0.30))
sapelo_bat_p1_boot

sapelo_bat_p1_boot<-as.ggplot(sapelo_bat_p1_boot)
sapelo_bat_p1_boot



#Sapelovirus p2
africa_sapelo_p2_bootscan <- read.csv(file = "africa_sapelo_p2_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_sapelo_p2_bootscan)

#move to long
long.sim_nt <- melt(africa_sapelo_p2_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","NC_033820",
                                                                                        "OQ818344","LC508226","OM104039", "OQ818342"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "E. dupreanum sapelovirus KEL272"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "R. madagascariensis sapelovirus MIZ243"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820"] <- "Eidolon helvum sapelovirus: NC_033820"
long.sim_nt$accession[long.sim_nt$accession == "OQ818344"] <- "E. dupreanum sapelovirus KEL298"
long.sim_nt$accession[long.sim_nt$accession == "OQ818342"] <- "E. dupreanum sapelovirus KEL273"
long.sim_nt$accession[long.sim_nt$accession == "LC508226"] <- "Porcine sapelovirus: LC508226"
long.sim_nt$accession[long.sim_nt$accession == "OM104039"] <- "Porcine sapelovirus: OM104039"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum sapelovirus KEL272",
                                                                  "R. madagascariensis sapelovirus MIZ243",
                                                                  "Eidolon helvum sapelovirus: NC_033820",
                                                                  "E. dupreanum sapelovirus KEL298",
                                                                  "E. dupreanum sapelovirus KEL273",
                                                                  "Porcine sapelovirus: LC508226",
                                                                  "Porcine sapelovirus: OM104039"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: E. dupreanum sapelovirus KEL233"))

sapelovirus_bat_p2_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
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
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_bat_p2_boot

#put gene map with PySimPlot
sapelo_bat_p2_boot<-sapelovirus_bat_p2_boot/ictv_sapelo_p2+plot_layout(nrow=2,  heights = c(2, 0.30))
sapelo_bat_p2_boot

sapelo_bat_p2_boot<-as.ggplot(sapelo_bat_p2_boot)
sapelo_bat_p2_boot



#Sapovirus all
africa_caliciviridae_all_bootscan <- read.csv(file = "africa_caliciviridae_all_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_caliciviridae_all_bootscan)

#move to long
long.sim_nt <- melt(africa_caliciviridae_all_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818340","OQ818345","OQ818347","OQ818348",
                                                                                    "NC_033776","OM105025"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818340"] <- "E. dupreanum sapovirus KEL207"
long.sim_nt$accession[long.sim_nt$accession == "OQ818345"] <- "R. madagascariensis sapovirus MIZ179"
long.sim_nt$accession[long.sim_nt$accession == "OQ818347"] <- "R. madagascariensis sapovirus MIZ219"
long.sim_nt$accession[long.sim_nt$accession == "OQ818348"] <- "R. madagascariensis sapovirus MIZ345"
long.sim_nt$accession[long.sim_nt$accession == "NC_033776"] <- "Eidolon helvum sapovirus: NC_033776"
long.sim_nt$accession[long.sim_nt$accession == "OM105025"] <- "Shellfish unclassified caliciviridae: OM105025"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum sapovirus KEL207",
                                                                  "R. madagascariensis sapovirus MIZ179",
                                                                  "R. madagascariensis sapovirus MIZ219",
                                                                  "R. madagascariensis sapovirus MIZ345",
                                                                  "Eidolon helvum sapovirus: NC_033776",
                                                                  "Shellfish unclassified caliciviridae: OM105025"))
#and plot
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: E. dupreanum sapovirus KEL166"))

sapovirus_all_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
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
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_all_bat_boot

#put gene map with PySimPlot
sapo_all_boot<-sapovirus_all_bat_boot/ictv_sapo_all+plot_layout(nrow=2,  heights = c(1, 0.30))
sapo_all_boot

sapo_all_boot<-as.ggplot(sapo_all_boot)
sapo_all_boot




#Sapovirus full
africa_caliciviridae_full_bootscan <- read.csv(file = "africa_caliciviridae_full_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_caliciviridae_full_bootscan)

#move to long
long.sim_nt <- melt(africa_caliciviridae_full_bootscan, id.vars = c("pointer"), measure.vars = c("JN699046",
                                                                                                "NC_033776","OM105025"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "JN699046"] <- "Human norovirus: JN699046"
long.sim_nt$accession[long.sim_nt$accession == "NC_033776"] <- "Eidolon helvum sapovirus: NC_033776"
long.sim_nt$accession[long.sim_nt$accession == "OM105025"] <- "Shellfish unclassified caliciviridae: OM105025"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Human norovirus: JN699046",
                                                                  "Eidolon helvum sapovirus: NC_033776",
                                                                  "Shellfish unclassified caliciviridae: OM105025"))
#and plot
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: E. dupreanum sapovirus KEL166"))

sapovirus_full_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
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
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_full_bat_boot

#put gene map with PySimPlot
sapo_full_boot<-sapovirus_full_bat_boot/ictv_sapo_full+plot_layout(nrow=2,  heights = c(1, 0.30))
sapo_full_boot

sapo_full_boot<-as.ggplot(sapo_full_boot)
sapo_full_boot



#Teschovirus
africa_tescho_bootscan <- read.csv(file = "africa_tescho_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_tescho_bootscan)

#move to long
long.sim_nt <- melt(africa_tescho_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818318","OQ818324","OM966657","OM105029"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818318"] <- "E. dupreanum teschovirus KEL164"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324"] <- "R. madagascariensis teschovirus MIZ205"
long.sim_nt$accession[long.sim_nt$accession == "OM966657"] <- "Sus scrofa teschovirus: OM966657"
long.sim_nt$accession[long.sim_nt$accession == "OM105029"] <- "Sus scrofa teschovirus: OM105029"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum teschovirus KEL164",
                                                                  "R. madagascariensis teschovirus MIZ205",
                                                                  "Sus scrofa teschovirus: OM966657",
                                                                  "Sus scrofa teschovirus: OM105029"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: R. madagascariensis teschovirus MIZ190"))

teschovirus_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
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
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

teschovirus_bat_boot

#put gene map with PySimPlot
tescho_boot<-teschovirus_bat_boot/ictv_tescho+plot_layout(nrow=2,  heights = c(1, 0.30))
tescho_boot

tescho_boot<-as.ggplot(tescho_boot)
tescho_boot


##Now put the whole figure together


#all bootscan together for supplementary: 
bootscan_all<-plot_grid(kun_bat_all_boot,
                        tescho_boot,
                           sapelo_bat_full_boot,
                    nrow=1,
                    labels="AUTO",  label_size = 23, align = "hv", axis="b")
bootscan_all

#export 20x5 inch PDF landscape

