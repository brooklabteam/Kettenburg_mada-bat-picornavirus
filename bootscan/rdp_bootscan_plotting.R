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
homewd="/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/bootscan/"
setwd("~/Desktop/developer/mada-bat-picornavirus/bootscan/gene_maps")

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
setwd("~/Desktop/developer/mada-bat-picornavirus/bootscan/rdp_bootscan_csv")

#colzpalette<-c("#F8766D","#C49A00","#53B400","#A58AFF","#00B6EB","darkorange1","#FB61D7")
colzpalette<-c("darkorange1","#00B6EB","#FB61D7")


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

long.sim_nt$accession[long.sim_nt$accession == "Bat_picornavirus_3...OQ818328"] <- "Bat picornavirus 3 (major parent) - OQ818328* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "Bat_picornavirus_3...PP766469"] <- "Bat picornavirus 3 (major parent) - PP766469* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818328...PP766469"] <- "OQ818328* (minor parent) - PP766469* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Bat picornavirus 3 (major parent) - OQ818328* (minor parent)",
                                                                  "Bat picornavirus 3 (major parent) - PP766469* (recombinant)",
                                                                  "OQ818328* (minor parent) - PP766469* (recombinant)"))

long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste("Reference: R. madagascariensis picornavirus 3: PP766469*"))

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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_map_nt

#only bootscan sig 1.339x10-2 from 2670-3533bp, bootstrap support vertical cutoff of 69
#Bat picornavirus 3 used to infer unknown parent

#put gene map with bootscan
bat_picorna_nt<-batpicorna_map_nt/map_batpicorna+plot_layout(nrow=2,  heights = c(1, 0.2))
bat_picorna_nt

bat_picorna_nt<-as.ggplot(bat_picorna_nt)
bat_picorna_nt

#Using PP766469 as reference seq 
batpicorna_map_left <- read.csv(file = "batpicorna_rdp_align_left_PP766469plot.csv", header = T, stringsAsFactors = F)
head(batpicorna_map_left)

#move to long
long.sim_nt <- melt(batpicorna_map_left, id.vars = c("pointer"), measure.vars = c("PP766471...OQ818325","PP766471...PP766469",
                                                                             "OQ818325...PP766469"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "PP766471...OQ818325"] <- "PP766471* (major parent) - OQ818325* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "PP766471...PP766469"] <- "PP766471* (major parent) - PP766469* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818325...PP766469"] <- "OQ818325* (minor parent) - PP766469* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("PP766471* (major parent) - OQ818325* (minor parent)",
                                                                  "PP766471* (major parent) - PP766469* (recombinant)",
                                                                  "OQ818325* (minor parent) - PP766469* (recombinant)"))

long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste("Reference: R. madagascariensis picornavirus 3: PP766469*"))

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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_map_nt_left

#bootscan sig 4.395x10-2 from 1375-1534, bootstrap support vertical cutoff of 69, chimaera sig 2.386x10-2

#put gene map with bootscan
bat_picorna_nt_left<-batpicorna_map_nt_left/map_batpicorna_left+plot_layout(nrow=2,  heights = c(1, 0.2))
bat_picorna_nt_left

bat_picorna_nt_left<-as.ggplot(bat_picorna_nt_left)
bat_picorna_nt_left

#Using shanbavirus A as reference seq 
batpicorna_map_mid <- read.csv(file = "batpicorna_rdp_align_mid_shanbavirusAplot.csv", header = T, stringsAsFactors = F)
head(batpicorna_map_mid)

#move to long
long.sim_nt <- melt(batpicorna_map_mid, id.vars = c("pointer"), measure.vars = c("PP766469...OQ818325","PP766469...Shanbavirus_A",
                                                                                  "OQ818325...Shanbavirus_A"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "PP766469...OQ818325"] <- "PP766469* (major parent) - OQ818325* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "PP766469...Shanbavirus_A"] <- "PP766469* (major parent) - Shanbavirus_A (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818325...Shanbavirus_A"] <- "OQ818325* (minor parent) - Shanbavirus_A (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("PP766469* (major parent) - OQ818325* (minor parent)",
                                                                  "PP766469* (major parent) - Shanbavirus_A (recombinant)",
                                                                  "OQ818325* (minor parent) - Shanbavirus_A (recombinant)"))

long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste("Reference: Shanbavirus A"))

#Plot nucleotide
batpicorna_map_nt_mid <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Bootstrap support")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "mid",
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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_map_nt_mid

#only bootscan sig 2.478x10-2 from 1540-2311bp, bootstrap support vertical cutoff of 69
#OQ818328 used to infer unknown parent

#put gene map with bootscan
bat_picorna_nt_mid<-batpicorna_map_nt_mid/map_batpicorna_mid+plot_layout(nrow=2,  heights = c(1, 0.2))
bat_picorna_nt_mid

bat_picorna_nt_mid<-as.ggplot(bat_picorna_nt_mid)
bat_picorna_nt_mid




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

long.sim_nt$accession[long.sim_nt$accession == "NC_028981.1...PP766455"] <- "NC_028981.1 (major parent) - PP766455* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "NC_028981.1...NC_028366.1"] <- "NC_028981.1 (major parent) - NC_028366.1 (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "PP766455...NC_028366.1"] <- "PP766455* (minor parent) - NC_028366.1 (recombinant)"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("NC_028981.1 (major parent) - PP766455* (minor parent)",
                                                                  "NC_028981.1 (major parent) - NC_028366.1 (recombinant)",
                                                                  "PP766455* (minor parent) - NC_028366.1 (recombinant)"))

long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: Hepatovirus H2"))

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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepato_map_nt

#bootscan sig 1.472x10-8 from 1-790 and 7740-7917bp, GENECONV sig 8.691x10-07, MaxChi sig 3.983x10-2, Chimaera sig 3.803x10-03
#nc_028981_1 used to infer unknown parent

#put gene map with bootscan
hepato_nt<-hepato_map_nt/map_hepato+plot_layout(nrow=2,  heights = c(1, 0.2))
hepato_nt

hepato_nt<-as.ggplot(hepato_nt)
hepato_nt




#Kobuvirus
#Using OP287812.1 as reference
kobuvirus_nt_map <- read.csv(file = "kobu_rdp_align_OP287812.1plot.csv", header = T, stringsAsFactors = F) 
head(kobuvirus_nt_map)

#move to long
long.sim_nt <- melt(kobuvirus_nt_map, id.vars = c("pointer"), measure.vars = c("Rhinolophus_kobuvirus...NC_034971.1","Rhinolophus_kobuvirus...OP287812.1",
                                                                               "NC_034971.1...OP287812.1"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "Rhinolophus_kobuvirus...NC_034971.1"] <- "Rhinolophus kobuvirus (major parent) - NC_034971.1 (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "Rhinolophus_kobuvirus...OP287812.1"] <- "Rhinolophus kobuvirus (major parent) - OP287812.1 (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "NC_034971.1...OP287812.1"] <- "NC_034971.1 (minor parent) - OP287812.1 (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rhinolophus kobuvirus (major parent) - NC_034971.1 (minor parent)",
                                                                  "Rhinolophus kobuvirus (major parent) - OP287812.1 (recombinant)",
                                                                  "NC_034971.1 (minor parent) - OP287812.1 (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. dupreanum kobuvirus: OP287812.1"))

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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobu_map_nt

#bootscan sig 1.064x10-3 from 4669 to 6403bp, maxchi sig 1.118x10-03, chimaera sig 4.388x10-03
#NC_034971_1 used to infer unknown parent, same recombination event seen in OQ818322 and PP766456

#put gene map with bootscan
kobu_nt<-kobu_map_nt/map_kobu+plot_layout(nrow=2,  heights = c(1, 0.2))
kobu_nt

kobu_nt<-as.ggplot(kobu_nt)
kobu_nt

#Using OP287812.1 as reference
kobuvirus_nt_map_mid <- read.csv(file = "kobu_rdp_align_mid_OP287812.1plot.csv", header = T, stringsAsFactors = F) 
head(kobuvirus_nt_map_mid)

#move to long
long.sim_nt <- melt(kobuvirus_nt_map_mid, id.vars = c("pointer"), measure.vars = c("NC_034971.1...Rhinolophus_kobuvirus",
                                                                                   "NC_034971.1...OP287812.1",
                                                                               "Rhinolophus_kobuvirus...OP287812.1"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_034971.1...Rhinolophus_kobuvirus"] <- "NC_034971.1 (major parent) - Rhinolophus kobuvirus (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "NC_034971.1...OP287812.1"] <- "NC_034971.1 (major parent) - OP287812.1 (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "Rhinolophus_kobuvirus...OP287812.1"] <- "Rhinolophus kobuvirus (minor parent) - OP287812.1 (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("NC_034971.1 (major parent) - Rhinolophus kobuvirus (minor parent)",
                                                                  "NC_034971.1 (major parent) - OP287812.1 (recombinant)",
                                                                  "Rhinolophus kobuvirus (minor parent) - OP287812.1 (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. dupreanum kobuvirus: OP287812.1"))

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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobu_map_nt_mid

#Bootscan sig 2.591x10-01 from 1419 to 2872, GENECONV sig 3.380x10-01, MaxChi sig 2.654x10-03, Chimaera sig 4.015x10-03
#Rhinolophus kpbuvirus used to infer unknown parent...all novel seqs have same recombination event

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

long.sim_nt$accession[long.sim_nt$accession == "Mischivirus_B...NC_026470.1"] <- "Mischivirus B (major parent) - NC_026470.1 (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "Mischivirus_B...OQ818316"] <- "Mischivirus B (major parent) - OQ818316* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "NC_026470.1...OQ818316"] <- "NC_026470.1 (minor parent) - OQ818316* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Mischivirus B (major parent) - NC_026470.1 (minor parent)",
                                                                  "Mischivirus B (major parent) - OQ818316* (recombinant)",
                                                                  "NC_026470.1 (minor parent) - OQ818316* (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: P. rufus mischivirus: OQ818316*"))

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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischi_map_nt

#put gene map with bootscan
mischi_nt<-mischi_map_nt/map_mischi+plot_layout(nrow=2,  heights = c(1, 0.2))
mischi_nt

mischi_nt<-as.ggplot(mischi_nt)
mischi_nt




#Sapelovirus
#reference sequence OQ818321
sapelovirus_full_nt_map <- read.csv(file = "sapelo_rdp_align_OQ818321plot.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapelovirus_full_nt_map)

#move to long
long.sim_nt <- melt(sapelovirus_full_nt_map, id.vars = c("pointer"), measure.vars = c("Pteropodidae_bat_sapelovirus...OQ818320","Pteropodidae_bat_sapelovirus...OQ818321",
                                                                                      "OQ818320...OQ818321"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "Pteropodidae_bat_sapelovirus...OQ818320"] <- "Pteropodidae bat sapelovirus (major parent) - OQ818320* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "Pteropodidae_bat_sapelovirus...OQ818321"] <- "Pteropodidae bat sapelovirus (major parent) - OQ818321* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818320...OQ818321"] <- "OQ818320* (minor parent) - OQ818321* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Pteropodidae bat sapelovirus (major parent) - OQ818320* (minor parent)",
                                                                  "Pteropodidae bat sapelovirus (major parent) - OQ818321* (recombinant)",
                                                                  "OQ818320* (minor parent) - OQ818321* (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: R. madagascariensis sapelovirus 2: OQ818321*"))

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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelo_map_nt

#put gene map with bootscan
sapelo_nt<-sapelo_map_nt/map_sapelo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapelo_nt

sapelo_nt<-as.ggplot(sapelo_nt)
sapelo_nt




#Teschovirus
#Reference of OQ818318
teschovirus_nt_map <- read.csv(file = "tescho_rdp_align_OQ818318plot.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_nt_map)

#move to long
long.sim_nt <- melt(teschovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818324...LC386158.1","OQ818324...OQ818318",
                                                                                 "LC386158.1...OQ818318"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818324...LC386158.1"] <- "OQ818324* (major parent) - LC386158.1 (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324...OQ818318"] <- "OQ818324* (major parent) - OQ818318* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "LC386158.1...OQ818318"] <- "LC386158.1 (minor parent) - OQ818318* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("OQ818324* (major parent) - LC386158.1 (minor parent)",
                                                                  "OQ818324* (major parent) - OQ818318* (recombinant)",
                                                                  "LC386158.1 (minor parent) - OQ818318* (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: E. dupreanum teschovirus 1: OQ818318*"))

tescho_map_nt1 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

tescho_map_nt1

#put gene map with bootscan
tescho_nt1<-tescho_map_nt1/map_tescho+plot_layout(nrow=2,  heights = c(1, 0.2))
tescho_nt1

tescho_nt1<-as.ggplot(tescho_nt1)
tescho_nt1

#reference sequence of Pteropodidae bat teschovirus
teschovirus_nt_map <- read.csv(file = "tescho_rdp_align_Pteropodidae_bat_teschovirus_plot.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_nt_map)

#move to long
long.sim_nt <- melt(teschovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818323...OQ818318","OQ818323...Pteropodidae_bat_teschovirus",
                                                                                 "OQ818318...Pteropodidae_bat_teschovirus"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818323...OQ818318"] <- "OQ818323* (major parent) - OQ818318* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818323...Pteropodidae_bat_teschovirus"] <- "OQ818323* (major parent) - Pteropodidae bat teschovirus (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818318...Pteropodidae_bat_teschovirus"] <- "OQ818318* (minor parent) - Pteropodidae bat teschovirus (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("OQ818323* (major parent) - OQ818318* (minor parent)",
                                                                  "OQ818323* (major parent) - Pteropodidae bat teschovirus (recombinant)",
                                                                  "OQ818318* (minor parent) - Pteropodidae bat teschovirus (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: Pteropodidae bat teschovirus"))

tescho_map_nt2 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

tescho_map_nt2

#put gene map with bootscan
tescho_nt2<-tescho_map_nt2/map_tescho+plot_layout(nrow=2,  heights = c(1, 0.2))
tescho_nt2

tescho_nt2<-as.ggplot(tescho_nt2)
tescho_nt2

#Reference of LC386158.1
teschovirus_nt_all_map <- read.csv(file = "tescho_rdp_align_all_LC386158.1plot.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_nt_all_map)

#move to long
long.sim_nt <- melt(teschovirus_nt_all_map, id.vars = c("pointer"), measure.vars = c("Pteropodidae_bat_teschovirus...OQ818324","Pteropodidae_bat_teschovirus...LC386158.1",
                                                                                 "OQ818324...LC386158.1"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "Pteropodidae_bat_teschovirus...OQ818324"] <- "Pteropodidae bat teschovirus (major parent) - OQ818324* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "Pteropodidae_bat_teschovirus...LC386158.1"] <- "Pteropodidae bat teschovirus (major parent) - LC386158.1 (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324...LC386158.1"] <- "OQ818324* (minor parent) - LC386158.1 (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Pteropodidae bat teschovirus (major parent) - OQ818324* (minor parent)",
                                                                  "Pteropodidae bat teschovirus (major parent) - LC386158.1 (recombinant)",
                                                                  "OQ818324* (minor parent) - LC386158.1 (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: Porcine teschovirus 16: LC386158.1"))

tescho_map_all_nt1 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

tescho_map_all_nt1

#put gene map with bootscan
tescho_all_nt1<-tescho_map_all_nt1/map_tescho_all+plot_layout(nrow=2,  heights = c(1, 0.2))
tescho_all_nt1

tescho_all_nt1<-as.ggplot(tescho_all_nt1)
tescho_all_nt1

#reference sequence of Pteropodidae bat teschovirus
teschovirus_nt_all_map <- read.csv(file = "tescho_rdp_align_Pteropodidae_bat_teschovirus_plot.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_nt_all_map)

#move to long
long.sim_nt <- melt(teschovirus_nt_all_map, id.vars = c("pointer"), measure.vars = c("OQ818323...OQ818318","OQ818323...Pteropodidae_bat_teschovirus",
                                                                                 "OQ818318...Pteropodidae_bat_teschovirus"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818323...OQ818318"] <- "OQ818323* (major parent) - OQ818318* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818323...Pteropodidae_bat_teschovirus"] <- "OQ818323* (major parent) - Pteropodidae bat teschovirus (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818318...Pteropodidae_bat_teschovirus"] <- "OQ818318* (minor parent) - Pteropodidae bat teschovirus (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("OQ818323* (major parent) - OQ818318* (minor parent)",
                                                                  "OQ818323* (major parent) - Pteropodidae bat teschovirus (recombinant)",
                                                                  "OQ818318* (minor parent) - Pteropodidae bat teschovirus (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: Pteropodidae bat teschovirus"))

tescho_map_all_nt2 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

tescho_map_all_nt2

#put gene map with bootscan
tescho_all_nt2<-tescho_map_all_nt2/map_tescho_all+plot_layout(nrow=2,  heights = c(1, 0.2))
tescho_all_nt2

tescho_all_nt2<-as.ggplot(tescho_all_nt2)
tescho_all_nt2




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

long.sim_nt$accession[long.sim_nt$accession == "R_leschenaulti_sapovirus...PP766459"] <- "R. leschenaulti sapovirus (major parent) - PP766459* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "R_leschenaulti_sapovirus...E_helvum_sapovirus"] <- "R. leschenaulti sapovirus (major parent) - E. helvum sapovirus (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "PP766459...E_helvum_sapovirus"] <- "PP766459* (minor parent) - E. helvum sapovirus (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. leschenaulti sapovirus (major parent) - PP766459* (minor parent)",
                                                                  "R. leschenaulti sapovirus (major parent) - E. helvum sapovirus (recombinant)",
                                                                  "PP766459* (minor parent) - E. helvum sapovirus (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. helvum sapovirus"))

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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapo_map_nt

#put gene map with bootscan
sapo_nt<-sapo_map_nt/map_sapo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapo_nt

sapo_nt<-as.ggplot(sapo_nt)
sapo_nt

#reference of PP766470
sapovirus_nt_map <- read.csv(file = "sapo_rdp_align_mid_PP766470plot.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapovirus_nt_map)

#move to long
long.sim_nt <- melt(sapovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818347...ON872527.1","OQ818347...PP766470",
                                                                               "ON872527.1...PP766470"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818347...ON872527.1"] <- "OQ818347* (major parent) - ON872527.1 (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818347...PP766470"] <- "OQ818347 *(major parent) - PP766470* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "ON872527.1...PP766470"] <- "ON872527.1 (minor parent) - PP766470* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("OQ818347* (major parent) - ON872527.1 (minor parent)",
                                                                  "OQ818347 *(major parent) - PP766470* (recombinant)",
                                                                  "ON872527.1 (minor parent) - PP766470* (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: R. madagascariensis sapovirus 2: PP766470*"))

sapo_map_nt_mid <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapo_map_nt_mid

#put gene map with bootscan
sapo_nt_mid<-sapo_map_nt_mid/map_sapo_mid+plot_layout(nrow=2,  heights = c(1, 0.2))
sapo_nt_mid

sapo_nt_mid<-as.ggplot(sapo_nt_mid)
sapo_nt_mid

#reference of E. helvum sapovirus
sapovirus_nt_map <- read.csv(file = "sapo_rdp_align_right_H_helvum_sapovirus_plot.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapovirus_nt_map)

#move to long
long.sim_nt <- melt(sapovirus_nt_map, id.vars = c("pointer"), measure.vars = c("PP766460...PP766477","PP766460...E_helvum_sapovirus",
                                                                               "PP766477...E_helvum_sapovirus"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "PP766460...PP766477"] <- "PP766460* (major parent) - PP766477* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "PP766460...E_helvum_sapovirus"] <- "PP766460* (major parent) - E. helvum sapovirus (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "PP766477...E_helvum_sapovirus"] <- "PP766477* (minor parent) - E. helvum sapovirus (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("PP766460* (major parent) - PP766477* (minor parent)",
                                                                  "PP766460* (major parent) - E. helvum sapovirus (recombinant)",
                                                                  "PP766477* (minor parent) - E. helvum sapovirus (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. helvum sapovirus"))

sapo_map_nt_right <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapo_map_nt_right

#put gene map with bootscan
sapo_nt_right<-sapo_map_nt_right/map_sapo_right+plot_layout(nrow=2,  heights = c(1, 0.2))
sapo_nt_right

sapo_nt_right<-as.ggplot(sapo_nt_right)
sapo_nt_right



##Now put the whole figure together
bootscan<-plot_grid(bat_picorna_nt,bat_picorna_nt_left, hepato_nt, kobu_nt, kobu_nt_mid,
                    mischi_nt, sapelo_nt, sapo_nt, sapo_nt_right, sapo_nt_right, 
                    tescho_nt1, tescho_nt2, tescho_nt_all_nt1, tescho_all_nt2,
                      ncol=4,
                      labels="AUTO", label_size = 23, align = "hv", axis="b")
bootscan


#export 20x23 landscape PDF

