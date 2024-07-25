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


#This is to make a figure of the representative bootscan figures from recombination detection program with their matching genome plots

##First make the genome plots
homewd="/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/recombination/"
setwd("~/Desktop/developer/mada-bat-picornavirus/recombination/gene_maps")

#Load the gene data
map <- read.csv("bootscan_alignment_genes.csv", header = T, stringsAsFactors = F)
map$gene<-factor(map$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                        "VP1", "2A", "2B", "2C", "3A", "3B",
                                        "3C", "3D", "NS1/NS2","Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Minor structural protein","3'UTR"))
#Load the peptide files
map_pep <- read.csv("bootscan_alignment_peptides.csv", header = T, stringsAsFactors = F)
map_pep$gene<-factor(map_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                              "VP1", "2A", "2B", "2C", "3A", "3B",
                                              "3C", "3D", "NS1/NS2","Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Minor structural protein","3'UTR"))
#Load the feature file in case its needed
map_feat <- read.csv("bootscan_alignment_features.csv", header = T, stringsAsFactors = F)
map_feat$gene<-factor(map_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                "VP1", "2A", "2B", "2C", "3A", "3B",
                                                "3C", "3D", "NS1/NS2","Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Minor structural protein","3'UTR"))

#Pick colors for genes
colz=c("5'UTR"="gold", "L"="royalblue","VP4"="paleturquoise3", "VP2"="skyblue1", "VP0"="royalblue4", "VP3"="steelblue1",
       "VP1"="cadetblue1", "2A"="orange1", "2B"="sienna2", "2C"="darkorange1", "3A"="palevioletred1", "3B"="plum",
       "3C"="rosybrown1", "3D"="pink3", "Helicase"="darkseagreen1","NS4"="darkolivegreen1","Vpg"="seagreen1","Pro-Pol"="palegreen1","NS1/NS2"="darkgreen",
       "Polyprotein"="azure3",
       "Minor structural protein"="black","3'UTR"="gold")

colz2=c("5'UTR"="grey", "L"="white","VP4"="white", "VP2"="white", "VP0"="white", "VP3"="white",
       "VP1"="white", "2A"="white", "2B"="white", "2C"="white", "3A"="white", "3B"="white",
       "3C"="white", "3D"="white", "Helicase"="white","NS4"="white","Vpg"="white","Pro-Pol"="white","NS1/NS2"="white",
       "Polyprotein"="black",
       "Minor structural protein"="black","3'UTR"="grey")


#Plot map and BLAST together plots

#Subset map data by virus
map_batpicornavirus<-subset(map,molecule=="Bat picornavirus")
map_batpicornavirus_left<-subset(map,molecule=="Bat picornavirus left")
map_batpicornavirus_mid<-subset(map,molecule=="Bat picornavirus mid")
map_batpicornavirus_right<-subset(map,molecule=="Bat picornavirus right")
map_hepatovirus<-subset(map,molecule=="Hepatovirus")
map_kobuvirus<-subset(map,molecule=="Kobuvirus")
map_kobuvirus_mid<-subset(map,molecule=="Kobuvirus mid")
map_mischivirus<-subset(map,molecule=="Mischivirus")
map_sapovirus<-subset(map,molecule=="Sapovirus")
map_sapovirus_mid<-subset(map,molecule=="Sapovirus mid")
map_sapovirus_right<-subset(map,molecule=="Sapovirus right")
map_sapelovirus<-subset(map,molecule=="Sapelovirus")
map_teschovirus<-subset(map,molecule=="Teschovirus")
map_teschovirus_all<-subset(map,molecule=="Teschovirus all")

map_batpicornavirus_pep<-subset(map_pep,molecule=="Bat picornavirus")
map_batpicornavirus_pep_left<-subset(map_pep,molecule=="Bat picornavirus left")
map_batpicornavirus_pep_mid<-subset(map_pep,molecule=="Bat picornavirus mid")
map_batpicornavirus_pep_right<-subset(map_pep,molecule=="Bat picornavirus right")
map_hepatovirus_pep<-subset(map_pep,molecule=="Hepatovirus")
map_kobuvirus_pep<-subset(map_pep,molecule=="Kobuvirus")
map_kobuvirus_pep_mid<-subset(map_pep,molecule=="Kobuvirus mid")
map_mischivirus_pep<-subset(map_pep,molecule=="Mischivirus")
map_sapovirus_pep<-subset(map_pep,molecule=="Sapovirus")
map_sapovirus_pep_mid<-subset(map_pep,molecule=="Sapovirus mid")
map_sapovirus_pep_right<-subset(map_pep,molecule=="Sapovirus right")
map_sapelovirus_pep<-subset(map_pep,molecule=="Sapelovirus")
map_teschovirus_pep<-subset(map_pep,molecule=="Teschovirus")
map_teschovirus_pep_all<-subset(map_pep,molecule=="Teschovirus all")

map_batpicornavirus_feat<-subset(map_feat,molecule=="Bat picornavirus")
map_batpicornavirus_feat_left<-subset(map_feat,molecule=="Bat picornavirus left")
map_batpicornavirus_feat_mid<-subset(map_feat,molecule=="Bat picornavirus mid")
map_batpicornavirus_feat_right<-subset(map_feat,molecule=="Bat picornavirus right")
map_hepatovirus_feat<-subset(map_feat,molecule=="Hepatovirus")
map_kobuvirus_feat<-subset(map_feat,molecule=="Kobuvirus")
map_kobuvirus_feat_mid<-subset(map_feat,molecule=="Kobuvirus mid")
map_mischivirus_feat<-subset(map_feat,molecule=="Mischivirus")
map_sapovirus_feat<-subset(map_feat,molecule=="Sapovirus")
map_sapovirus_feat_mid<-subset(map_feat,molecule=="Sapovirus mid")
map_sapovirus_feat_right<-subset(map_feat,molecule=="Sapovirus right")
map_sapelovirus_feat<-subset(map_feat,molecule=="Sapelovirus")
map_teschovirus_feat<-subset(map_feat,molecule=="Teschovirus")
map_teschovirus_feat_all<-subset(map_feat,molecule=="Teschovirus all")

#gene maps
map_batpicorna<-ggplot(map_batpicornavirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_batpicornavirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                  xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_batpicornavirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(0,8552),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_batpicorna



map_batpicorna_left<-ggplot(map_batpicornavirus_left, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_batpicornavirus_pep_left, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_batpicornavirus_feat_left,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(0,2835),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_batpicorna_left



map_batpicorna_mid<-ggplot(map_batpicornavirus_mid, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_batpicornavirus_pep_mid, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_batpicornavirus_feat_mid,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(2779,5425),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_batpicorna_mid


map_batpicorna_right<-ggplot(map_batpicornavirus_right, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_batpicornavirus_pep_right, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                     xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_batpicornavirus_feat_right,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(5528,6170),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_batpicorna_right


map_hepato<-ggplot(map_hepatovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_hepatovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_hepatovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(0,7920),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_hepato


map_kobu<-ggplot(map_kobuvirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_kobuvirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_kobuvirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_kobuvirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                            xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_kobuvirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(0,8535),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_kobu


map_kobu_mid<-ggplot(map_kobuvirus_mid, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_kobuvirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_kobuvirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_kobuvirus_pep_mid, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                           xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_kobuvirus_feat_mid,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(2485,5327),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_kobu_mid


map_mischi<-ggplot(map_mischivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_mischivirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_mischivirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_mischivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_mischivirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(0,8890),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_mischi


map_sapelo<-ggplot(map_sapelovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_sapelovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_sapelovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_sapelovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                   xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_sapelovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(0,8230),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_sapelo


map_tescho<-ggplot(map_teschovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_teschovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_teschovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_teschovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_teschovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(0,7700),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_tescho


map_tescho_all<-ggplot(map_teschovirus_all, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_teschovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_teschovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_teschovirus_pep_all, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                             xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_teschovirus_feat_all,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(4561,7593),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_tescho_all



map_sapo<-ggplot(map_sapovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_sapovirus_full_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_sapovirus_full_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_sapovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_sapovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(0,7920),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_sapo


map_sapo_mid<-ggplot(map_sapovirus_mid, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_sapovirus_full_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_sapovirus_full_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_sapovirus_pep_mid, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                           xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_sapovirus_feat_mid,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(2692,4850),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_sapo_mid


map_sapo_right<-ggplot(map_sapovirus_right, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_sapovirus_full_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_sapovirus_full_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_sapovirus_pep_right, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                           xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_sapovirus_feat_right,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(6666,7895),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_sapo_right


##Now get the plots for the bootscan
setwd("~/Desktop/developer/mada-bat-picornavirus/recombination/rdp_bootscan_csv_to_plot")

#colzpalette<-c("#F8766D","#C49A00","#53B400","#A58AFF","#00B6EB","darkorange1","#FB61D7")
colzpalette<-c("#3B9AB2","#EBCC2A","#F21A00")


##Bat picornavirus

#Using PP766469 as reference seq 
batpicorna_map <- read.csv(file = "batpicorna_rdp_align_PP766469plot.csv", header = T, stringsAsFactors = F)
head(batpicorna_map)

#move to long
long.sim_nt <- melt(batpicorna_map, id.vars = c("pointer"), measure.vars = c("Bat_picornavirus_3...OQ818328","Bat_picornavirus_3...PP766469",
                                                                             "OQ818328...PP766469"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "Bat_picornavirus_3...OQ818328"] <- "Bat picornavirus 3 (major parent) - R. madagascariensis picornavirus 1: OQ818328* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "Bat_picornavirus_3...PP766469"] <- "Bat picornavirus 3 (major parent) - R. madagascariensis picornavirus 3: PP766469* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818328...PP766469"] <- "R. madagascariensis picornavirus 1: OQ818328* (minor parent) - R. madagascariensis picornavirus 3: PP766469* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Bat picornavirus 3 (major parent) - R. madagascariensis picornavirus 1: OQ818328* (minor parent)",
                                                                  "Bat picornavirus 3 (major parent) - R. madagascariensis picornavirus 3: PP766469* (recombinant)",
                                                                  "R. madagascariensis picornavirus 1: OQ818328* (minor parent) - R. madagascariensis picornavirus 3: PP766469* (recombinant)"))

long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste("Potential recombinant - R. madagascariensis picornavirus 3: PP766469*"))

#Plot nucleotide
batpicorna_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Bootstrap support")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  geom_hline(yintercept=0.69, linetype='dashed', col = 'black')+
  annotate("rect", xmin=528, xmax=688, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_map_nt

#RDP sig 4.164x10-02, 3Seq sig 2.150x10-02 from 528 to 688bp
#OP963617.1 used to infer unknown parent

#put gene map with bootscan
bat_picorna_nt<-batpicorna_map_nt/map_batpicorna+plot_layout(nrow=2,  heights = c(1, 0.2))
bat_picorna_nt

bat_picorna_nt<-as.ggplot(bat_picorna_nt)
bat_picorna_nt

#Using PP766469 as reference seq 
batpicorna_map_left <- read.csv(file = "batpicorna_rdp_align_left_PP766469plot.csv", header = T, stringsAsFactors = F)
head(batpicorna_map_left)

#move to long
long.sim_nt <- melt(batpicorna_map_left, id.vars = c("pointer"), measure.vars = c("OQ818328...Bat_picornavirus_3","OQ818328...PP766469",
                                                                             "Bat_picornavirus_3...PP766469"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818328...Bat_picornavirus_3"] <- "R. madagascariensis picornavirus 1: OQ818328* (major parent) - Bat picornavirus 3* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818328...PP766469"] <- "R. madagascariensis picornavirus 1: OQ818328* (major parent) - R. madagascariensis picornavirus 3: PP766469* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "Bat_picornavirus_3...PP766469"] <- "Bat picornavirus 3 (minor parent) - R. madagascariensis picornavirus 3: PP766469* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. madagascariensis picornavirus 1: OQ818328* (major parent) - Bat picornavirus 3* (minor parent)",
                                                                  "R. madagascariensis picornavirus 1: OQ818328* (major parent) - R. madagascariensis picornavirus 3: PP766469* (recombinant)",
                                                                  "Bat picornavirus 3 (minor parent) - R. madagascariensis picornavirus 3: PP766469* (recombinant)"))

long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste("Potential recombinant - R. madagascariensis picornavirus 3: PP766469*"))

#Plot nucleotide
batpicorna_map_nt_left <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Bootstrap support")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  geom_hline(yintercept=0.69, linetype='dashed', col = 'black')+
  annotate("rect", xmin=102, xmax=264, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_map_nt_left

#RDP sig 2.392x10-04, bootscan sig 1.25x10-04, 3Seq sig 3.48x10-02 from 102-264bp

#put gene map with bootscan
bat_picorna_nt_left<-batpicorna_map_nt_left/map_batpicorna_left+plot_layout(nrow=2,  heights = c(1, 0.2))
bat_picorna_nt_left

bat_picorna_nt_left<-as.ggplot(bat_picorna_nt_left)
bat_picorna_nt_left

#Using OQ818346 as reference seq 
batpicorna_map_right <- read.csv(file = "batpicorna_rdp_align_right_OQ818346plot.csv", header = T, stringsAsFactors = F)
head(batpicorna_map_right)

#move to long
long.sim_nt <- melt(batpicorna_map_right, id.vars = c("pointer"), measure.vars = c("OP963617.1...Bat_picornavirus_3","OP963617.1...OQ818346",
                                                                                  "Bat_picornavirus_3...OQ818346"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OP963617.1...Bat_picornavirus_3"] <- "Unknown (inferred by Bat picornavirus BtSY4) (major parent) - Bat picornavirus 3 (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "OP963617.1...OQ818346"] <- "Unknown (inferred by Bat picornavirus BtSY4) (major parent) - R. madagascariensis picornavirus 2: OQ818346* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "Bat_picornavirus_3...OQ818346"] <- "Bat picornavirus 3 (minor parent) - R. madagascariensis picornavirus 2: OQ818346* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Unknown (inferred by Bat picornavirus BtSY4) (major parent) - Bat picornavirus 3 (minor parent)",
                                                                  "Unknown (inferred by Bat picornavirus BtSY4) (major parent) - R. madagascariensis picornavirus 2: OQ818346* (recombinant)",
                                                                  "Bat picornavirus 3 (minor parent) - R. madagascariensis picornavirus 2: OQ818346* (recombinant)"))

long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste("Potential recombinant - R. madagascariensis picornavirus 2: OQ818346*"))

#Plot nucleotide
batpicorna_map_nt_right <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Bootstrap support")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) +
  geom_hline(yintercept=0.69, linetype='dashed', col = 'black')+
  annotate("rect", xmin=1, xmax=100, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=516, xmax=635, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,200,400,600),
                     labels = c(0,5730,5930,6130),expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_map_nt_right

#RDP sig 1.307x10-01, MaxChi sig 1.201x10-03, Chimaera sif 1.143x10-02, 3Seq sig 3.686x10-02 from 516-635nt, 1-100nt region 2
#OP963617.1 used to infer unknown parent

#put gene map with bootscan
bat_picorna_nt_right<-batpicorna_map_nt_right/map_batpicorna_right+plot_layout(nrow=2,  heights = c(1, 0.2))
bat_picorna_nt_right

bat_picorna_nt_right<-as.ggplot(bat_picorna_nt_right)
bat_picorna_nt_right




#Hepatovirus
# Using NC_028366.1 as reference
hepato_map_nt <- read.csv(file = "hepato_rdp_align_NC_028366.1plot.csv", header = T, stringsAsFactors = F) #animo acid
head(hepato_map_nt)

#move to long
long.sim_nt <- melt(hepato_map_nt, id.vars = c("pointer"), measure.vars = c("NC_028981.1...PP766455","NC_028981.1...NC_028366.1",
                                                                            "PP766455...NC_028366.1"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_028981.1...PP766455"] <- "Unknown (inferred by Tupaia hepatovirus) (major parent) - E. dupreanum hepatovirus: PP766455* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "NC_028981.1...NC_028366.1"] <- "Unknown (inferred by Tupaia hepatovirus)  (major parent) - E. helvum hepatovirus M32Eidhel2010 (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "PP766455...NC_028366.1"] <- "E. dupreanum hepatovirus: PP766455* (minor parent) - E. helvum hepatovirus M32Eidhel2010 (recombinant)"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Unknown (inferred by Tupaia hepatovirus) (major parent) - E. dupreanum hepatovirus: PP766455* (minor parent)",
                                                                  "Unknown (inferred by Tupaia hepatovirus)  (major parent) - E. helvum hepatovirus M32Eidhel2010 (recombinant)",
                                                                  "E. dupreanum hepatovirus: PP766455* (minor parent) - E. helvum hepatovirus M32Eidhel2010 (recombinant)"))

long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Potential recombinant - E. helvum hepatovirus M32Eidhel2010"))

hepato_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Bootstrap support")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  geom_hline(yintercept=0.69, linetype='dashed', col = 'black')+
  annotate("rect", xmin=1, xmax=790, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=7740, xmax=7917, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=3495, xmax=3885, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepato_map_nt

#RDP sig 1.201x10-04, Geneconv sig 8.691x10-07, Bootscan sig 1.472x10-08,   MaxChi sig 3.983x10-02, Chimaera sig 3.803x10-03 from 1-790 then 7740 to 7917bp
#maybe sig also from 3495 to 3885bp
#NC_028981.1 used to infer unknown parent

#put gene map with bootscan
hepato_nt<-hepato_map_nt/map_hepato+plot_layout(nrow=2,  heights = c(1, 0.2))
hepato_nt

hepato_nt<-as.ggplot(hepato_nt)
hepato_nt




#Kobuvirus
#Using OQ818322 as reference
kobuvirus_nt_map <- read.csv(file = "kobu_rdp_align_OQ818322plot.csv", header = T, stringsAsFactors = F) 
head(kobuvirus_nt_map)

#move to long
long.sim_nt <- melt(kobuvirus_nt_map, id.vars = c("pointer"), measure.vars = c("Rhinolophus_kobuvirus...NC_034971.1","Rhinolophus_kobuvirus...OQ818322",
                                                                               "NC_034971.1...OQ818322"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "Rhinolophus_kobuvirus...NC_034971.1"] <- "Rhinolophus kobuvirus (major parent) - Unknown (inferred by canine kobuvirus) (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "Rhinolophus_kobuvirus...OQ818322"] <- "Rhinolophus kobuvirus (major parent) - E. dupreanum kobuvirus: OQ818322* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "NC_034971.1...OQ818322"] <- "Unknown (inferred by canine kobuvirus) (minor parent) - E. dupreanum kobuvirus: OQ818322* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rhinolophus kobuvirus (major parent) - Unknown (inferred by canine kobuvirus) (minor parent)",
                                                                  "Rhinolophus kobuvirus (major parent) - E. dupreanum kobuvirus: OQ818322* (recombinant)",
                                                                  "Unknown (inferred by canine kobuvirus) (minor parent) - E. dupreanum kobuvirus: OQ818322* (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Potential recombinant: E. dupreanum kobuvirus: OQ818322*"))

kobu_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Bootstrap support")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  geom_hline(yintercept=0.69, linetype='dashed', col = 'black')+
  annotate("rect", xmin=4669, xmax=5714, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=3338, xmax=3695, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobu_map_nt

#RDP sig 5.097x10-03, bootscan sig 1.294x10-01, MaxChi sig 4.368x10-03, Chimaera sig 7.185x10-03 from 4669 to 5714bp for region 2
#RDP sig 9.155x10-05, bootscan sig 1.533x10-02, MaxChi sig 3.61x10-02, Chimaera sig 8.473x10-03, 3Seq sig 4.378x10-04 from 3338 to 3695bp
#nc_034971.1 used to infer unknown parent, OQ818322, PP766456, and OP287812.1 have the same recombination pattern

#put gene map with bootscan
kobu_nt<-kobu_map_nt/map_kobu+plot_layout(nrow=2,  heights = c(1, 0.2))
kobu_nt

kobu_nt<-as.ggplot(kobu_nt)
kobu_nt

#Using OQ818322 as reference
kobuvirus_nt_map_mid <- read.csv(file = "kobu_rdp_align_mid_OQ818322plot.csv", header = T, stringsAsFactors = F) 
head(kobuvirus_nt_map_mid)

#move to long
long.sim_nt <- melt(kobuvirus_nt_map_mid, id.vars = c("pointer"), measure.vars = c("NC_034971.1...Rhinolophus_kobuvirus",
                                                                                   "NC_034971.1...OQ818322",
                                                                               "Rhinolophus_kobuvirus...OQ818322"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_034971.1...Rhinolophus_kobuvirus"] <- "Canine kobuvirus (major parent) - E. dupreanum kobuvirus: OQ818322* (inferred by Rhinolophus kobuvirus) (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "NC_034971.1...OQ818322"] <- "Canine kobuvirus (major parent) -  E. dupreanum kobuvirus: OQ818322* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "Rhinolophus_kobuvirus...OQ818322"] <- "Unknown (inferred by Rhinolophus kobuvirus) (minor parent) - E. dupreanum kobuvirus: OQ818322* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Canine kobuvirus (major parent) - E. dupreanum kobuvirus: OQ818322* (inferred by Rhinolophus kobuvirus) (minor parent)",
                                                                  "Canine kobuvirus (major parent) -  E. dupreanum kobuvirus: OQ818322* (recombinant)",
                                                                  "Unknown (inferred by Rhinolophus kobuvirus) (minor parent) - E. dupreanum kobuvirus: OQ818322* (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Potential recombinant - E. dupreanum kobuvirus: OQ818322*"))

kobu_map_nt_mid <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Bootstrap support")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  geom_hline(yintercept=0.69, linetype='dashed', col = 'black')+
  annotate("rect", xmin=2208, xmax=2872, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,1000,2000),
                     labels = c(0,3489,4489),expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobu_map_nt_mid

#RDP sig 3.435x10-02, Bootscan sig 2.591x10-01, MaxChi sig 2.654x10-03, Chimaera sig 4.015x10-03 from 2208-2872bp
#Rhinolophus kobuvirus used to infer unknown parent, same recombination event in all Malagasy seq

#put gene map with bootscan
kobu_nt_mid<-kobu_map_nt_mid/map_kobu_mid+plot_layout(nrow=2,  heights = c(1, 0.2))
kobu_nt_mid

kobu_nt_mid<-as.ggplot(kobu_nt_mid)
kobu_nt_mid




#Mischivirus
#Using OQ818316 as reference
mischivirus_nt_map <- read.csv(file = "mischi_rdp_align_OQ818316plot.csv", header = T, stringsAsFactors = F) #Amino acid
head(mischivirus_nt_map)

#move to long
long.sim_nt <- melt(mischivirus_nt_map, id.vars = c("pointer"), measure.vars = c("Mischivirus_B...NC_026470.1","Mischivirus_B...OQ818316",
                                                                                 "NC_026470.1...OQ818316"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "Mischivirus_B...NC_026470.1"] <- "Unknown (inferred by Mischivirus B) (major parent) - African bat icavirus A (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "Mischivirus_B...OQ818316"] <- "Unknown (inferred by Mischivirus B) (major parent) - P. rufus mischivirus: OQ818316* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "NC_026470.1...OQ818316"] <- "African bat icavirus A (minor parent) - P. rufus mischivirus: OQ818316* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Unknown (inferred by Mischivirus B) (major parent) - African bat icavirus A (minor parent)",
                                                                  "Unknown (inferred by Mischivirus B) (major parent) - P. rufus mischivirus: OQ818316* (recombinant)",
                                                                  "African bat icavirus A (minor parent) - P. rufus mischivirus: OQ818316* (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Potential recombinant - P. rufus mischivirus: OQ818316*"))

mischi_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Bootstrap support")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  geom_hline(yintercept=0.69, linetype='dashed', col = 'black')+
  annotate("rect", xmin=6151, xmax=6526, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischi_map_nt

#RDP sig 1.643x10-02, MaxChi sig 3.020x10-02 from 6151 to 6526bp
#Mischivirus B used to infer unknown parent

#put gene map with bootscan
mischi_nt<-mischi_map_nt/map_mischi+plot_layout(nrow=2,  heights = c(1, 0.2))
mischi_nt

mischi_nt<-as.ggplot(mischi_nt)
mischi_nt




#Sapelovirus
#reference sequence OQ818329
sapelovirus_full_nt_map <- read.csv(file = "sapelo_rdp_align_OQ818329plot.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapelovirus_full_nt_map)

#move to long
long.sim_nt <- melt(sapelovirus_full_nt_map, id.vars = c("pointer"), measure.vars = c("NC_003987.1...NC_033820.1","NC_003987.1...OQ818329",
                                                                                      "NC_033820.1...OQ818329"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_003987.1...NC_033820.1"] <- "Unknown (inferred by porcine sapelovirus 1) (major parent) - Bat sapelovirus Bat/CAM/Sap-P24/2013 (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "NC_003987.1...OQ818329"] <- "Unknown (inferred by porcine sapelovirus 1) (major parent) - R. madagascariensis sapelovirus 1: OQ818329* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820.1...OQ818329"] <- "Bat sapelovirus Bat/CAM/Sap-P24/2013 (minor parent) - R. madagascariensis sapelovirus 1: OQ818329* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Unknown (inferred by porcine sapelovirus 1) (major parent) - Bat sapelovirus Bat/CAM/Sap-P24/2013 (minor parent)",
                                                                  "Unknown (inferred by porcine sapelovirus 1) (major parent) - R. madagascariensis sapelovirus 1: OQ818329* (recombinant)",
                                                                  "Bat sapelovirus Bat/CAM/Sap-P24/2013 (minor parent) - R. madagascariensis sapelovirus 1: OQ818329* (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Potential recombinant - R. madagascariensis sapelovirus 1: OQ818329*"))

sapelo_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Bootstrap support")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  geom_hline(yintercept=0.69, linetype='dashed', col = 'black')+
  annotate("rect", xmin=159, xmax=554, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=4170, xmax=5092, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelo_map_nt

#RDP sig 7.434x10-02, GENECONV sig 7.391X10-03, Bootscan sig 1.409x10-03, MaxChi sig 2.652x10-04, Chimaera sig 2.641x10-02 from 159-554 in region 1
#RDP sig 1.606x10-03, GENECONV sig 3.194x10-04, MaxChi sig 2.189x10-02, Chimaera sig 2.713x10-02 from 4710 to 5092
#NC_003987.1 used to infer unknown parent

#put gene map with bootscan
sapelo_nt<-sapelo_map_nt/map_sapelo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapelo_nt

sapelo_nt<-as.ggplot(sapelo_nt)
sapelo_nt




#Sapovirus
#reference of E. helvum sapovirus
sapovirus_nt_map <- read.csv(file = "sapo_rdp_align_E_helvum_sapovirus_plot.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapovirus_nt_map)

#move to long
long.sim_nt <- melt(sapovirus_nt_map, id.vars = c("pointer"), measure.vars = c("R_leschenaulti_sapovirus...PP766459","R_leschenaulti_sapovirus...E_helvum_sapovirus",
                                                                               "PP766459...E_helvum_sapovirus"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "R_leschenaulti_sapovirus...PP766459"] <- "Unknown (inferred by R. leschenaulti sapovirus) (major parent) - E. dupreanum sapovirus 1: PP766459* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "R_leschenaulti_sapovirus...E_helvum_sapovirus"] <- "Unknown (inferred by R. leschenaulti sapovirus) (major parent) - E. helvum sapovirus (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "PP766459...E_helvum_sapovirus"] <- "E. dupreanum sapovirus 1: PP766459* (minor parent) - E. helvum sapovirus (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Unknown (inferred by R. leschenaulti sapovirus) (major parent) - E. dupreanum sapovirus 1: PP766459* (minor parent)",
                                                                  "Unknown (inferred by R. leschenaulti sapovirus) (major parent) - E. helvum sapovirus (recombinant)",
                                                                  "E. dupreanum sapovirus 1: PP766459* (minor parent) - E. helvum sapovirus (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Potential recombinant - E. helvum sapovirus"))

sapo_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Bootstrap support")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_blank(),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  geom_hline(yintercept=0.69, linetype='dashed', col = 'black')+
  annotate("rect", xmin=4100, xmax=5872, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapo_map_nt

#RDP sig 1.234x10-05, GENECONV sig 1.620x10-07, Bootscan sig 3.122x10-05, MaxChi sig 2.976x10-02, Chimaera sig 7.021x10-03, 3Seq sig 3.006x10-04 from 4100-5872bp
#RDP sig 6.218x10-04 from 1800-2623
#R. leschenaulti used to infer unknown parent

#put gene map with bootscan
sapo_nt<-sapo_map_nt/map_sapo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapo_nt

sapo_nt<-as.ggplot(sapo_nt)
sapo_nt



#To plot figure 4
fig4<-plot_grid(kobu_nt,kobu_nt_mid, sapo_nt,
                ncol=3,
                labels="AUTO", label_size = 23, align = "hv", axis="b")
fig4
#only showing kobuvirus and sapovirus since they have the strongest support for recombination


#To plot supplementary for everything else
recombination_supp<-plot_grid(mischi_nt, hepato_nt,bat_picorna_nt, bat_picorna_nt_left,
                           bat_picorna_nt_right, sapelo_nt,
                           ncol=3,
                           labels="AUTO", label_size = 23, align = "hv", axis="b")
recombination_supp




# save figs

homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 

ggsave(file = paste0(homewd, "/final_figures/Fig4_recombination.pdf"),
       plot = fig4,
       units="mm",  
       width=150, 
       height=40, 
       scale=4, 
       dpi=300)

ggsave(file = paste0(homewd, "/final_figures/supplemental/SfigX_recombination.pdf"),
       plot = recombination_supp,
       units="mm",  
       width=160, 
       height=80, 
       scale=4, 
       dpi=300)


