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
map_hepatovirus<-subset(map,molecule=="Hepatovirus")
map_kobuvirus<-subset(map,molecule=="Kobuvirus")
map_kunsagivirus<-subset(map,molecule=="Kunsagivirus")
map_mischivirus<-subset(map,molecule=="Mischivirus")
map_sapovirus<-subset(map,molecule=="Sapovirus")
map_sapelovirus<-subset(map,molecule=="Sapelovirus")
map_teschovirus<-subset(map,molecule=="Teschovirus")

map_batpicornavirus_pep<-subset(map_pep,molecule=="Bat picornavirus")
map_hepatovirus_pep<-subset(map_pep,molecule=="Hepatovirus")
map_kobuvirus_pep<-subset(map_pep,molecule=="Kobuvirus")
map_kunsagivirus_pep<-subset(map_pep,molecule=="Kunsagivirus")
map_mischivirus_pep<-subset(map_pep,molecule=="Mischivirus")
map_sapovirus_pep<-subset(map_pep,molecule=="Sapovirus")
map_sapelovirus_pep<-subset(map_pep,molecule=="Sapelovirus")
map_teschovirus_pep<-subset(map_pep,molecule=="Teschovirus")

map_batpicornavirus_feat<-subset(map_feat,molecule=="Bat picornavirus")
map_hepatovirus_feat<-subset(map_feat,molecule=="Hepatovirus")
map_kobuvirus_feat<-subset(map_feat,molecule=="Kobuvirus")
map_kunsagivirus_feat<-subset(map_feat,molecule=="Kunsagivirus")
map_mischivirus_feat<-subset(map_feat,molecule=="Mischivirus")
map_sapovirus_feat<-subset(map_feat,molecule=="Sapovirus")
map_sapelovirus_feat<-subset(map_feat,molecule=="Sapelovirus")
map_teschovirus_feat<-subset(map_feat,molecule=="Teschovirus")

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


map_kun<-ggplot(map_kunsagivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_kunsagivirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_kunsagivirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_kunsagivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                               xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_kunsagivirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(0,7975),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_kun


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



##Now get the plots for the bootscan
setwd("~/Desktop/developer/mada-bat-picornavirus/bootscan/rdp_bootscan_csv")

#colzpalette<-c("#F8766D","#C49A00","#53B400","#A58AFF","#00B6EB","darkorange1","#FB61D7")
colzpalette<-c("darkorange1","#00B6EB","#FB61D7")


##Bat picornavirus

##from reference
batpicorna_map <- read.csv(file = "batpicorna_bootscan_align_nt.csv", header = T, stringsAsFactors = F)
head(batpicorna_map)

#move to long
long.sim_nt <- melt(batpicorna_map, id.vars = c("pointer"), measure.vars = c("OQ818328","PP766469","Bat_picornavirus_3",
                                                                             "Shanbavirus_A","OP963617.1"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818328"] <- "R. madagascariensis picornavirus 1: OQ818328*"
long.sim_nt$accession[long.sim_nt$accession == "PP766469"] <- "R. madagascariensis picornavirus 3: PP766469*"
long.sim_nt$accession[long.sim_nt$accession == "Bat_picornavirus_3"] <- "Bat picornavirus 3"
long.sim_nt$accession[long.sim_nt$accession == "Shanbavirus_A"] <- "Shanbavirus A"
long.sim_nt$accession[long.sim_nt$accession == "OP963617.1"] <- "Bat picornavirus BtSY4"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. madagascariensis picornavirus 1: OQ818328*",
                                                                  "R. madagascariensis picornavirus 3: PP766469*",
                                                                  "Bat picornavirus 3",
                                                                  "Shanbavirus A",
                                                                  "Bat picornavirus BtSY4"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste("Reference: R. madagascariensis picornavirus 1: OQ818325*"))

#Plot nucleotide
batpicorna_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Permuted trees")+xlab("Genome position")+
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

#put gene map with bootscan
bat_picorna_nt<-batpicorna_map_nt/map_batpicorna+plot_layout(nrow=2,  heights = c(1, 0.2))
bat_picorna_nt

bat_picorna_nt<-as.ggplot(bat_picorna_nt)
bat_picorna_nt




#Hepatovirus
#nucleotide
hepato_map_nt <- read.csv(file = "hepato_bootscan_align_nt.csv", header = T, stringsAsFactors = F) #animo acid
head(hepato_map_nt)

#move to long
long.sim_nt <- melt(hepato_map_nt, id.vars = c("pointer"), measure.vars = c("PP766457","Hipposideros_bat_hepatovirus","NC_028366.1",
                                                                            "NC_028981.1","NC_038313.1"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "PP766457"] <- "E. dupreanum hepatovirus: PP766457*"
long.sim_nt$accession[long.sim_nt$accession == "Hipposideros_bat_hepatovirus"] <- "Hipposideros bat hepatovirus"
long.sim_nt$accession[long.sim_nt$accession == "NC_028366.1"] <- "Hepatovirus H2"
long.sim_nt$accession[long.sim_nt$accession == "NC_028981.1"] <- "Tupaia hepatovirus A"
long.sim_nt$accession[long.sim_nt$accession == "NC_038313.1"] <- "Bat hepatovirus SMG18520Minmav2014"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum hepatovirus: PP766457*",
                                                                  "Hipposideros bat hepatovirus",
                                                                  "Hepatovirus H2",
                                                                  "Tupaia hepatovirus A",
                                                                  "Bat hepatovirus SMG18520Minmav2014"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. dupreanum hepatovirus: PP766455*"))

hepato_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Permuted trees")+xlab("Genome position")+
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

#put gene map with bootscan
hepato_nt<-hepato_map_nt/map_hepato+plot_layout(nrow=2,  heights = c(1, 0.2))
hepato_nt

hepato_nt<-as.ggplot(hepato_nt)
hepato_nt




#Kobuvirus
#Nucleotide
kobuvirus_nt_map <- read.csv(file = "kobu_bootscan_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(kobuvirus_nt_map)

#move to long
long.sim_nt <- melt(kobuvirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818322","PP766456","NC_034971.1",
                                                                               "Rhinolophus_kobuvirus"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818322"] <- "E. dupreanum kobuvirus: OQ818322*"
long.sim_nt$accession[long.sim_nt$accession == "PP766456"] <- "E. dupreanum kobuvirus: PP766456*"
long.sim_nt$accession[long.sim_nt$accession == "NC_034971.1"] <- "Canine kobuvirus"
long.sim_nt$accession[long.sim_nt$accession == "Rhinolophus_kobuvirus"] <- "Rhinolophus kobuvirus"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum kobuvirus: OQ818322*",
                                                                  "E. dupreanum kobuvirus: PP766456*",
                                                                  "Canine kobuvirus",
                                                                  "Rhinolophus kobuvirus"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. dupreanum kobuvirus: OP287812.1"))

kobu_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Permuted trees")+xlab("Genome position")+
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
  guides(colour = guide_legend(nrow = 2))+
  scale_color_manual(values=colzpalette) + 
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobu_map_nt

#put gene map with bootscan
kobu_nt<-kobu_map_nt/map_kobu+plot_layout(nrow=2,  heights = c(1, 0.2))
kobu_nt

kobu_nt<-as.ggplot(kobu_nt)
kobu_nt




#kunsagivirus
#Nucleotide
kunsagivirus_nt_map <- read.csv(file = "kunsagi_bootscan_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(kunsagivirus_nt_map)

#move to long
long.sim_nt <- melt(kunsagivirus_nt_map, id.vars = c("pointer"), measure.vars = c("NC_033818.1","OP589993.1","NC_034206.1"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_033818.1"] <- "Kunsagivirus B"
long.sim_nt$accession[long.sim_nt$accession == "OP589993.1"] <- "Auskunsag virus"
long.sim_nt$accession[long.sim_nt$accession == "NC_034206.1"] <- "Bakunsa virus"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Kunsagivirus B", "Auskunsag virus","Bakunsa virus"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. dupreanum kunsagivirus: OQ818317*"))

kunsagi_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Permuted trees")+xlab("Genome position")+
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
  guides(colour = guide_legend(nrow = 1))+
  scale_color_manual(values=colzpalette) + 
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagi_map_nt

#put gene map with bootscan
kunsagi_nt<-kunsagi_map_nt/map_kun+plot_layout(nrow=2,  heights = c(1, 0.2))
kunsagi_nt

kunsagi_nt<-as.ggplot(kunsagi_nt)
kunsagi_nt




#Mischivirus
#Plot nucleotide
mischivirus_nt_map <- read.csv(file = "mischi_bootscan_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(mischivirus_nt_map)

#move to long
long.sim_nt <- melt(mischivirus_nt_map, id.vars = c("pointer"), measure.vars = c("Bat_mischivirus_5","Mischivirus_B",
                                                                                 "NC_026470.1"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "Bat_mischivirus_5"] <- "Bat mischivirus 5"
long.sim_nt$accession[long.sim_nt$accession == "Mischivirus_B"] <- "Mischivirus B"
long.sim_nt$accession[long.sim_nt$accession == "NC_026470.1"] <- "Mischivirus A"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Bat mischivirus 5", "Mischivirus B",
                                                                  "Mischivirus A"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: P. rufus mischivirus: OQ818316*"))

mischi_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Permuted trees")+xlab("Genome position")+
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
  guides(colour = guide_legend(nrow = 1))+
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
#Plot nucleotide
sapelovirus_full_nt_map <- read.csv(file = "sapelo_bootscan_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapelovirus_full_nt_map)

#move to long
long.sim_nt <- melt(sapelovirus_full_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818320","NC_003987.1","NC_033820.1",
                                                                                      "Pteropodidae_bat_sapelovirus"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "E. dupreanum sapelovirus 2: OQ818321*"
long.sim_nt$accession[long.sim_nt$accession == "OQ818320"] <- "E. dupreanum sapelovirus 1: OQ818320*"
long.sim_nt$accession[long.sim_nt$accession == "NC_003987.1"] <- "Porcine sapelovirus 1"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820.1"] <- "Bat sapelovirus"
long.sim_nt$accession[long.sim_nt$accession == "Pteropodidae_bat_sapelovirus"] <- "Pteropodidae bat sapelovirus"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum sapelovirus 1: OQ818320*",
                                                                  "E. dupreanum sapelovirus 2: OQ818321*",
                                                                  "Porcine sapelovirus 1",
                                                                  "Bat sapelovirus",
                                                                  "Pteropodidae bat sapelovirus"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: R. madagascariensis sapelovirus 1: OQ818329*"))

sapelo_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Permuted trees")+xlab("Genome position")+
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




#Teschovirus 1
#Plot nucleotide
teschovirus_nt_map <- read.csv(file = "tescho_bootscan_align_nt1.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_nt_map)

#move to long
long.sim_nt <- melt(teschovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818318","OQ818323","LC386158.1",
                                                                                 "Pteropodidae_bat_teschovirus"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818318"] <- "E. dupreanum teschovirus 1: OQ818318*"
long.sim_nt$accession[long.sim_nt$accession == "OQ818323"] <- "R. madagascariensis teschovirus 1: OQ818323*"
long.sim_nt$accession[long.sim_nt$accession == "LC386158.1"] <- "Porcine teschovirus 16"
long.sim_nt$accession[long.sim_nt$accession == "Pteropodidae_bat_teschovirus"] <- "Pteropodidae bat teschovirus"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum teschovirus 1: OQ818318*",
                                                                  "R. madagascariensis teschovirus 1: OQ818323*",
                                                                  "Porcine teschovirus 16",
                                                                  "Pteropodidae bat teschovirus"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: R. madagascariensis teschovirus 2: OQ818324*"))

tescho_map_nt1 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Permuted trees")+xlab("Genome position")+
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
  guides(colour = guide_legend(nrow = 2))+
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




#Teschovirus 2
#Plot nucleotide
teschovirus_nt_map <- read.csv(file = "tescho_bootscan_align_nt2.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_nt_map)

#move to long
long.sim_nt <- melt(teschovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818323","OQ818324","LC386158.1",
                                                                                 "Pteropodidae_bat_teschovirus"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818324"] <- "R. madagascariensis teschovirus 2: OQ818324*"
long.sim_nt$accession[long.sim_nt$accession == "OQ818323"] <- "R. madagascariensis teschovirus 1: OQ818323*"
long.sim_nt$accession[long.sim_nt$accession == "LC386158.1"] <- "Porcine teschovirus 16"
long.sim_nt$accession[long.sim_nt$accession == "Pteropodidae_bat_teschovirus"] <- "Pteropodidae bat teschovirus"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. madagascariensis teschovirus 1: OQ818323*",
                                                                  "R. madagascariensis teschovirus 2: OQ818324*",
                                                                  "Porcine teschovirus 16",
                                                                  "Pteropodidae bat teschovirus"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: E. dupreanum teschovirus 1: OQ818318*"))

tescho_map_nt2 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Permuted trees")+xlab("Genome position")+
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
  guides(colour = guide_legend(nrow = 2))+
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




#Sapovirus
#Plot nucleotide
sapovirus_nt_map <- read.csv(file = "sapo_bootscan_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapovirus_nt_map)

#move to long
long.sim_nt <- melt(sapovirus_nt_map, id.vars = c("pointer"), measure.vars = c("E_helvum_sapovirus","OQ818319","R_leschenaulti_sapovirus",
                                                                               "ON872527.1"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "E_helvum_sapovirus"] <- "E. helvum sapovirus"
long.sim_nt$accession[long.sim_nt$accession == "OQ818319"] <- "E. dupreanum sapovirus 1: OQ818319*"
long.sim_nt$accession[long.sim_nt$accession == "R_leschenaulti_sapovirus"] <- "R. leschenaulti sapovirus"
long.sim_nt$accession[long.sim_nt$accession == "ON872527.1"] <- "Bat faecal sapovirus"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum sapovirus 1: OQ818319*",
                                                                  "E. helvum sapovirus",
                                                                  "R. leschenaulti sapovirus",
                                                                  "Bat faecal sapovirus"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. dupreanum sapovirus 1: PP766459*"))

sapo_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Permuted trees")+xlab("Genome position")+
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
  guides(colour = guide_legend(nrow = 2))+
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


##Now put the whole figure together
bootscan<-plot_grid(bat_picorna_nt,hepato_nt,kobu_nt,kunsagi_nt,
                      mischi_nt,sapelo_nt,sapo_nt, tescho_nt1, tescho_nt2,
                      ncol=3,
                      labels="AUTO", label_size = 23, align = "hv", axis="b")
bootscan


#export 20x23 landscape PDF

