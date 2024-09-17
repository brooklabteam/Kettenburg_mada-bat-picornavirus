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
homewd="/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/PySimPlot/"
setwd("~/Desktop/developer/mada-bat-picornavirus/PySimPlot/gene_maps")

#Load the gene data
map <- read.csv("pysimplot_alignment_genes.csv", header = T, stringsAsFactors = F)
map$gene<-factor(map$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                        "VP1","VP1 ", "2A", "2B", "2C", "3A", "3B",
                                        "3C", "3D", "NS1/NS2","Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Minor structural protein","3'UTR"))
#Load the peptide files
map_pep <- read.csv("pysimplot_alignment_peptides.csv", header = T, stringsAsFactors = F)
map_pep$gene<-factor(map_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                              "VP1", "VP1 ","2A", "2B", "2C", "3A", "3B",
                                              "3C", "3D", "NS1/NS2","Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Minor structural protein","3'UTR"))
#Load the feature file in case its needed
map_feat <- read.csv("pysimplot_alignment_features.csv", header = T, stringsAsFactors = F)
map_feat$gene<-factor(map_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                "VP1","VP1 ", "2A", "2B", "2C", "3A", "3B",
                                                "3C", "3D", "NS1/NS2","Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Minor structural protein","3'UTR"))

#Pick colors for genes
colz=c("5'UTR"="gold", "L"="royalblue","VP4"="paleturquoise3", "VP2"="skyblue1", "VP0"="royalblue4", "VP3"="steelblue1",
       "VP1"="cadetblue1", "2A"="orange1", "2B"="sienna2", "2C"="darkorange1", "3A"="palevioletred1", "3B"="plum",
       "3C"="rosybrown1", "3D"="pink3", "Helicase"="darkseagreen1","NS4"="darkolivegreen1","Vpg"="seagreen1","Pro-Pol"="palegreen1","NS1/NS2"="darkgreen",
       "Polyprotein"="azure3",
       "Minor structural protein"="black","3'UTR"="gold")

colz2=c("5'UTR"="cornflowerblue", "L"="white","VP4"="white", "VP2"="white", "VP0"="white", "VP3"="white",
       "VP1"="white", "VP1 "="orange1", "2A"="orange1", "2B"="white", "2C"="white", "3A"="orange1", "3B"="white",
       "3C"="white", "3D"="white", "Helicase"="white","NS4"="orange1","Vpg"="white","Pro-Pol"="white","NS1/NS2"="white",
       "Polyprotein"="black",
       "Minor structural protein"="black","3'UTR"="cornflowerblue")


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
map_cardiovirus<-subset(map,molecule=="Cardiovirus")

map_batpicornavirus_pep<-subset(map_pep,molecule=="Bat picornavirus")
map_hepatovirus_pep<-subset(map_pep,molecule=="Hepatovirus")
map_kobuvirus_pep<-subset(map_pep,molecule=="Kobuvirus")
map_kunsagivirus_pep<-subset(map_pep,molecule=="Kunsagivirus")
map_mischivirus_pep<-subset(map_pep,molecule=="Mischivirus")
map_sapovirus_pep<-subset(map_pep,molecule=="Sapovirus")
map_sapelovirus_pep<-subset(map_pep,molecule=="Sapelovirus")
map_teschovirus_pep<-subset(map_pep,molecule=="Teschovirus")
map_cardiovirus_pep<-subset(map_pep,molecule=="Cardiovirus")

map_batpicornavirus_feat<-subset(map_feat,molecule=="Bat picornavirus")
map_hepatovirus_feat<-subset(map_feat,molecule=="Hepatovirus")
map_kobuvirus_feat<-subset(map_feat,molecule=="Kobuvirus")
map_kunsagivirus_feat<-subset(map_feat,molecule=="Kunsagivirus")
map_mischivirus_feat<-subset(map_feat,molecule=="Mischivirus")
map_sapovirus_feat<-subset(map_feat,molecule=="Sapovirus")
map_sapelovirus_feat<-subset(map_feat,molecule=="Sapelovirus")
map_teschovirus_feat<-subset(map_feat,molecule=="Teschovirus")
map_cardiovirus_feat<-subset(map_feat,molecule=="Cardiovirus")

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
  scale_x_continuous(limits=c(0,8290),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_batpicorna


map_cardio<-ggplot(map_cardiovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=map_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=map_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = map_cardiovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black",
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=map_cardiovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz2)+
  theme_genes()+
  scale_x_continuous(limits=c(4000,7640),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_cardio


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
  scale_x_continuous(limits=c(0,7890),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,8420),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,7580),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,9080),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,7940),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,7445),expand=c(0,0))+
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
  scale_x_continuous(limits=c(0,7800),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
map_sapo



##Now get the plots for the PySimPlot, just the map_BLAST ones, there will be a separate PDF file with the table of African bat picorna similarities
setwd("~/Desktop/developer/mada-bat-picornavirus/PySimPlot/output")

#colzpalette<-c("#F8766D","#C49A00","#53B400","#A58AFF","#00B6EB","darkorange1","#FB61D7")
colzpalette<-c("#3B9AB2","#EBCC2A","#F21A00","darkslategray2","darkorange1","firebrick4")


##Bat picornavirus

##nt first
batpicorna_map <- read.csv(file = "batpicorna_align_nt.csv", header = T, stringsAsFactors = F)
head(batpicorna_map)

#move to long
long.sim_nt <- melt(batpicorna_map, id.vars = c("pointer"), measure.vars = c("PP766469","OQ818325","PP711912","PP711945","PP711930"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "PP766469"] <- "R. madagascariensis picornavirus 3: PP766469*"
long.sim_nt$accession[long.sim_nt$accession == "OQ818325"] <- "R. madagascariensis picornavirus 1: OQ818325*"
long.sim_nt$accession[long.sim_nt$accession == "PP711912"] <- "Rousettus bat picornavirus: PP711912"
long.sim_nt$accession[long.sim_nt$accession == "PP711945"] <- "Rousettus bat picornavirus: PP711945"
long.sim_nt$accession[long.sim_nt$accession == "PP711930"] <- "Rousettus bat picornavirus: PP711930"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. madagascariensis picornavirus 3: PP766469*", 
                                                                  "R. madagascariensis picornavirus 1: OQ818325*",
                                                                  "Rousettus bat picornavirus: PP711912",
                                                                  "Rousettus bat picornavirus: PP711945",
                                                                  "Rousettus bat picornavirus: PP711930"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

title<-expression(paste("Reference: R. madagascariensis picornavirus 1: OQ818325*"))

#Plot nucleotide
batpicorna_map_nt <- ggplot(long.sim_nt) + 
  annotate("rect", xmin=3168, xmax=3904, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5161, xmax=5579, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Nucleotide identity")+xlab("Genome position")+
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
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_map_nt

#put gene map with PySimPlot
bat_picorna_nt<-batpicorna_map_nt/map_batpicorna+plot_layout(nrow=2,  heights = c(1, 0.2))
bat_picorna_nt

bat_picorna_nt<-as.ggplot(bat_picorna_nt)
bat_picorna_nt

#then aa
batpicorna_map <- read.csv(file = "batpicorna_align_aa.csv", header = T, stringsAsFactors = F)
head(batpicorna_map)

#move to long
long.sim_aa <- melt(batpicorna_map, id.vars = c("pointer"), measure.vars = c("PP766469","OQ818325","PP711912","PP711945","PP711930"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "PP766469"] <- "R. madagascariensis picornavirus 3: PP766469*"
long.sim_aa$accession[long.sim_aa$accession == "OQ818325"] <- "R. madagascariensis picornavirus 1: OQ818325*"
long.sim_aa$accession[long.sim_aa$accession == "PP711912"] <- "Rousettus bat picornavirus: PP711912"
long.sim_aa$accession[long.sim_aa$accession == "PP711945"] <- "Rousettus bat picornavirus: PP711945"
long.sim_aa$accession[long.sim_aa$accession == "PP711930"] <- "Rousettus bat picornavirus: PP711930"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("R. madagascariensis picornavirus 3: PP766469*", 
                                                                  "R. madagascariensis picornavirus 1: OQ818325*",
                                                                  "Rousettus bat picornavirus: PP711912",
                                                                  "Rousettus bat picornavirus: PP711945",
                                                                  "Rousettus bat picornavirus: PP711930"))
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

title<-expression(paste("Reference: R. madagascariensis picornavirus 1: OQ818325*"))
## amino acid
batpicorna_map_aa <- ggplot(long.sim_aa) + 
  annotate("rect", xmin=3368/3.055, xmax=4104/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5261/3.055, xmax=5679/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Amino acid identity")+xlab("Genome position")+
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
  ggtitle(title)+
  scale_color_manual(values=colzpalette) + 
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_map_aa

#put gene map with PySimPlot
bat_picorna_aa<-batpicorna_map_aa/map_batpicorna+plot_layout(nrow=2,  heights = c(1, 0.2))
bat_picorna_aa

bat_picorna_aa<-as.ggplot(bat_picorna_aa)
bat_picorna_aa




#Hepatovirus
#nucleotide
hepato_map_nt <- read.csv(file = "hepato_align_nt.csv", header = T, stringsAsFactors = F) #animo acid
head(hepato_map_nt)

#move to long
long.sim_nt <- melt(hepato_map_nt, id.vars = c("pointer"), measure.vars = c("PP766457","KT452742","NC_028366","NC_028981"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "PP766457"] <- "E. dupreanum hepatovirus: PP766457*"
long.sim_nt$accession[long.sim_nt$accession == "KT452742"] <- "Bat hepatovirus: KT452742"
long.sim_nt$accession[long.sim_nt$accession == "NC_028366"] <- "Hepatovirus H2: NC_028366"
long.sim_nt$accession[long.sim_nt$accession == "NC_028981"] <- "Tupaia hepatovirus A: NC_028981"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum hepatovirus: PP766457*","Bat hepatovirus: KT452742",
                                                                  "Hepatovirus H2: NC_028366","Tupaia hepatovirus A: NC_028981"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100
#plot nucleotide
title<-expression(paste("Reference: E. dupreanum hepatovirus: PP766455*"))

hepato_map_nt <- ggplot(long.sim_nt) + 
  annotate("rect", xmin=3068, xmax=3804, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5061, xmax=5379, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Nucleotide identity")+xlab("Genome position")+
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
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepato_map_nt

#put gene map with PySimPlot
hepato_nt<-hepato_map_nt/map_hepato+plot_layout(nrow=2,  heights = c(1, 0.2))
hepato_nt

hepato_nt<-as.ggplot(hepato_nt)
hepato_nt

##amino acid
hepato_map_aa <- read.csv(file = "hepato_align_aa.csv", header = T, stringsAsFactors = F) #animo acid
head(hepato_map_aa)

#move to long
long.sim_aa <- melt(hepato_map_aa, id.vars = c("pointer"), measure.vars = c("PP766457","KT452742","NC_028366","NC_028981"))
unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "PP766457"] <- "E. dupreanum hepatovirus: PP766457*"
long.sim_aa$accession[long.sim_aa$accession == "KT452742"] <- "Bat hepatovirus: KT452742"
long.sim_aa$accession[long.sim_aa$accession == "NC_028366"] <- "Hepatovirus H2: NC_028366"
long.sim_aa$accession[long.sim_aa$accession == "NC_028981"] <- "Tupaia hepatovirus A: NC_028981"


long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("E. dupreanum hepatovirus: PP766457*","Bat hepatovirus: KT452742",
                                                                  "Hepatovirus H2: NC_028366","Tupaia hepatovirus A: NC_028981"))
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## Amino acid
title<-expression(paste("Reference: E. dupreanum hepatovirus: PP766455*"))

hepatovirus_map_aa <- ggplot(long.sim_aa) +
  annotate("rect", xmin=3068/3.055, xmax=3804/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5061/3.055, xmax=5379/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Amino acid identity")+xlab("Genome position")+
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
  #scale_fill_distiller()+
  ggtitle(title)+
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_map_aa

#put gene map with PySimPlot
hep_map_aa<-hepatovirus_map_aa/map_hepato+plot_layout(nrow=2,  heights = c(1, 0.2))
hep_map_aa

hep_map_aa<-as.ggplot(hep_map_aa)
hep_map_aa




#Kobuvirus
#Nucleotide
kobuvirus_nt_map <- read.csv(file = "kobu_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(kobuvirus_nt_map)

#move to long
long.sim_nt <- melt(kobuvirus_nt_map, id.vars = c("pointer"), measure.vars = c("PP766456","OP287812"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "PP766456"] <- "E. dupreanum kobuvirus: PP766456*"
long.sim_nt$accession[long.sim_nt$accession == "OP287812"] <- "E. dupreanum kobuvirus: OP287812"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum kobuvirus: PP766456*", "E. dupreanum kobuvirus: OP287812"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. dupreanum kobuvirus: OQ818322*"))

kobu_map_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Nucleotide identity")+xlab("Genome position")+
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
  annotate("rect", xmin=3768, xmax=4304, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5661, xmax=5979, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobu_map_nt

#put gene map with PySimPlot
kobu_nt<-kobu_map_nt/map_kobu+plot_layout(nrow=2,  heights = c(1, 0.2))
kobu_nt

kobu_nt<-as.ggplot(kobu_nt)
kobu_nt

#amino acid
kobuvirus_aa_map <- read.csv(file = "kobu_align_aa.csv", header = T, stringsAsFactors = F) #Amino acid
head(kobuvirus_aa_map)

#move to long
long.sim_aa <- melt(kobuvirus_aa_map, id.vars = c("pointer"), measure.vars = c("PP766456","OP287812"))
unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "PP766456"] <- "E. dupreanum kobuvirus: PP766456*"
long.sim_aa$accession[long.sim_aa$accession == "OP287812"] <- "E. dupreanum kobuvirus: OP287812"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("E. dupreanum kobuvirus: PP766456*", "E. dupreanum kobuvirus: OP287812"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## Amino acid
title<-expression(paste("Reference: E. dupreanum kobuvirus: OQ818322*"))

kobuvirus_map_aa <- ggplot(long.sim_aa) + 
  annotate("rect", xmin=3768/3.055, xmax=4304/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5661/3.055, xmax=5979/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Amino acid identity")+xlab("Genome position")+
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
  #scale_fill_distiller()+
  ggtitle(title)+
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000,4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobuvirus_map_aa

#put gene map with PySimPlot
kobu_map_aa<-kobuvirus_map_aa/map_kobu+plot_layout(nrow=2,  heights = c(1, 0.2))
kobu_map_aa

kobu_map_aa<-as.ggplot(kobu_map_aa)
kobu_map_aa




#kunsagivirus
#Nucleotide
kunsagivirus_nt_map <- read.csv(file = "kunsagi_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(kunsagivirus_nt_map)

#move to long
long.sim_nt <- melt(kunsagivirus_nt_map, id.vars = c("pointer"), measure.vars = c("NC_033818","OP589993"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_033818"] <- "Kunsagivirus B: NC_033818"
long.sim_nt$accession[long.sim_nt$accession == "OP589993"] <- "Auskunsag virus: OP589993"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Kunsagivirus B: NC_033818", "Auskunsag virus: OP589993"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. dupreanum kunsagivirus: OQ818317*"))

kunsagi_map_nt <- ggplot(long.sim_nt) + 
  annotate("rect", xmin=2768, xmax=3304, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=4861, xmax=5179, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Nucleotide identity")+xlab("Genome position")+
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
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagi_map_nt

#put gene map with PySimPlot
kunsagi_nt<-kunsagi_map_nt/map_kun+plot_layout(nrow=2,  heights = c(1, 0.2))
kunsagi_nt

kunsagi_nt<-as.ggplot(kunsagi_nt)
kunsagi_nt

#amino acid
kunsagivirus_aa_map <- read.csv(file = "kunsagi_align_aa.csv", header = T, stringsAsFactors = F) #Amino acid
head(kunsagivirus_aa_map)

#move to long
long.sim_aa <- melt(kunsagivirus_aa_map, id.vars = c("pointer"), measure.vars = c("NC_033818","OP589993"))
unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "NC_033818"] <- "Kunsagivirus B: NC_033818"
long.sim_aa$accession[long.sim_aa$accession == "OP589993"] <- "Auskunsag virus: OP589993"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Kunsagivirus B: NC_033818", "Auskunsag virus: OP589993"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## Amino acid
title<-expression(paste("Reference: E. dupreanum kunsagivirus: OQ818317*"))

kunsagivirus_map_aa <- ggplot(long.sim_aa) + 
  annotate("rect", xmin=2768/3.055, xmax=3304/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=4861/3.055, xmax=5179/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Amino acid identity")+xlab("Genome position")+
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
  #scale_fill_distiller()+
  ggtitle(title)+
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000,4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagivirus_map_aa

#put gene map with PySimPlot
kun_map_aa<-kunsagivirus_map_aa/map_kun+plot_layout(nrow=2,  heights = c(1, 0.2))
kun_map_aa

kun_map_aa<-as.ggplot(kun_map_aa)
kun_map_aa




#Mischivirus
#Plot nucleotide
mischivirus_nt_map <- read.csv(file = "mischi_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(mischivirus_nt_map)

#move to long
long.sim_nt <- melt(mischivirus_nt_map, id.vars = c("pointer"), measure.vars = c("JQ814851","MG888045",
                                                                                 "NC_026470"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "JQ814851"] <- "M. schreibersii picornavirus 1: JQ814851"
long.sim_nt$accession[long.sim_nt$accession == "MG888045"] <- "Mischivirus B: MG888045"
long.sim_nt$accession[long.sim_nt$accession == "NC_026470"] <- "Mischivirus A: NC_026470"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("M. schreibersii picornavirus 1: JQ814851", "Mischivirus B: MG888045",
                                                                  "Mischivirus A: NC_026470"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: P. rufus mischivirus: OQ818316*"))

mischi_map_nt <- ggplot(long.sim_nt) + 
  annotate("rect", xmin=4368, xmax=4704, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=6261, xmax=6779, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Nucleotide identity")+xlab("Genome position")+
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
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischi_map_nt

#put gene map with PySimPlot
mischi_nt<-mischi_map_nt/map_mischi+plot_layout(nrow=2,  heights = c(1, 0.2))
mischi_nt

mischi_nt<-as.ggplot(mischi_nt)
mischi_nt

#amino acid
mischivirus_aa_map <- read.csv(file = "mischi_align_aa.csv", header = T, stringsAsFactors = F) #Amino acid
head(mischivirus_aa_map)

#move to long
long.sim_aa <- melt(mischivirus_aa_map, id.vars = c("pointer"), measure.vars = c("JQ814851","MG888045",
                                                                                  "NC_026470"))
unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "JQ814851"] <- "M. schreibersii picornavirus 1: JQ814851"
long.sim_aa$accession[long.sim_aa$accession == "MG888045"] <- "Mischivirus B: MG888045"
long.sim_aa$accession[long.sim_aa$accession == "NC_026470"] <- "Mischivirus A: NC_026470"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("M. schreibersii picornavirus 1: JQ814851", "Mischivirus B: MG888045",
                                                                  "Mischivirus A: NC_026470"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## Amino acid
title<-expression(paste("Reference: P. rufus mischivirus: OQ818316*"))

mischivirus_map_aa <- ggplot(long.sim_aa) + 
  annotate("rect", xmin=4368/3.055, xmax=4704/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=6261/3.055, xmax=6779/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Amino acid identity")+xlab("Genome position")+
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
  #scale_fill_distiller()+
  ggtitle(title)+
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_map_aa

#put gene map with PySimPlot
mischi_map_aa<-mischivirus_map_aa/map_mischi+plot_layout(nrow=2,  heights = c(1, 0.2))
mischi_map_aa

mischi_map_aa<-as.ggplot(mischi_map_aa)
mischi_map_aa




#Sapelovirus
#Plot nucleotide
sapelovirus_full_nt_map <- read.csv(file = "sapelo_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapelovirus_full_nt_map)

#move to long
long.sim_nt <- melt(sapelovirus_full_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818320","OQ818329","PP711921","PP711911","NC_033820"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818320"] <- "E. dupreanum sapelovirus: OQ818320*"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "R. madagascariensis sapelovirus: OQ818329*"
long.sim_nt$accession[long.sim_nt$accession == "PP711921"] <- "Eidolon bat picornavirus: PP711921"
long.sim_nt$accession[long.sim_nt$accession == "PP711911"] <- "Rousettus bat picornavirus: PP711911"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820"] <- "Eidolon helvum sapelovirus: NC_033820"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum sapelovirus: OQ818320*",
                                                                  "R. madagascariensis sapelovirus: OQ818329*",
                                                                  "Rousettus bat picornavirus: PP711911",
                                                                  "Eidolon bat picornavirus: PP711921",
                                                                  "Eidolon helvum sapelovirus: NC_033820"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: E. dupreanum sapelovirus: OQ818321*"))

sapelo_map_nt <- ggplot(long.sim_nt) + 
  annotate("rect", xmin=3368, xmax=4004, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5461, xmax=5779, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Nucleotide identity")+xlab("Genome position")+
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
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelo_map_nt

#put gene map with PySimPlot
sapelo_nt<-sapelo_map_nt/map_sapelo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapelo_nt

sapelo_nt<-as.ggplot(sapelo_nt)
sapelo_nt

#amino acid
sapelovirus_full_aa_map <- read.csv(file = "sapelo_align_aa.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapelovirus_full_aa_map)

#move to long
long.sim_aa <- melt(sapelovirus_full_aa_map, id.vars = c("pointer"), measure.vars = c("OQ818320","OQ818329","PP711921","PP711911","NC_033820"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "OQ818320"] <- "E. dupreanum sapelovirus: OQ818320*"
long.sim_aa$accession[long.sim_aa$accession == "OQ818329"] <- "R. madagascariensis sapelovirus: OQ818329*"
long.sim_aa$accession[long.sim_aa$accession == "PP711921"] <- "Eidolon bat picornavirus: PP711921"
long.sim_aa$accession[long.sim_aa$accession == "PP711911"] <- "Rousettus bat picornavirus: PP711911"
long.sim_aa$accession[long.sim_aa$accession == "NC_033820"] <- "Eidolon helvum sapelovirus: NC_033820"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("E. dupreanum sapelovirus: OQ818320*",
                                                                  "R. madagascariensis sapelovirus: OQ818329*",
                                                                  "Rousettus bat picornavirus: PP711911",
                                                                  "Eidolon bat picornavirus: PP711921",
                                                                  "Eidolon helvum sapelovirus: NC_033820"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## nucleotide
title<-expression(paste("Reference: E. dupreanum sapelovirus: OQ818321*"))

sapelovirus_full_map_aa <- ggplot(long.sim_aa) + 
  annotate("rect", xmin=3368/3.055, xmax=4204/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5561/3.055, xmax=5879/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Amino acid identity")+xlab("Genome position")+
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
  #scale_fill_distiller()+
  ggtitle(title)+
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_full_map_aa

#put gene map with PySimPlot
sapelo_map_aa<-sapelovirus_full_map_aa/map_sapelo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapelo_map_aa

sapelo_map_aa<-as.ggplot(sapelo_map_aa)
sapelo_map_aa




#Teschovirus
#Plot nucleotide
teschovirus_nt_map <- read.csv(file = "tescho_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_nt_map)

#move to long
long.sim_nt <- melt(teschovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818324","OQ818318","PP711948",
                                                                                 "OR951334"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818324"] <- "R. madagascariensis teschovirus 2: OQ818324*"
long.sim_nt$accession[long.sim_nt$accession == "OQ818318"] <- "E. dupreanum teschovirus 1: OQ818318*"
long.sim_nt$accession[long.sim_nt$accession == "OR951334"] <- "Pteropodidae bat teschovirus: OR951334"
long.sim_nt$accession[long.sim_nt$accession == "PP711948"] <- "Rousettus bat picornavirus: PP711948"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. madagascariensis teschovirus 2: OQ818324*",
                                                                  "E. dupreanum teschovirus 1: OQ818318*",
                                                                  "Pteropodidae bat teschovirus: OR951334",
                                                                  "Rousettus bat picornavirus: PP711948"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Reference: R. madagascariensis teschovirus 1: OQ818323*"))
tescho_map_nt <- ggplot(long.sim_nt) + 
  annotate("rect", xmin=3368, xmax=3604, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=4961, xmax=5279, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Nucleotide identity")+xlab("Genome position")+
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
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

tescho_map_nt

#put gene map with PySimPlot
tescho_nt<-tescho_map_nt/map_tescho+plot_layout(nrow=2,  heights = c(1, 0.2))
tescho_nt

tescho_nt<-as.ggplot(tescho_nt)
tescho_nt

#amino acid
teschovirus_aa_map <- read.csv(file = "tescho_align_aa.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_aa_map)

#move to long
long.sim_aa <- melt(teschovirus_aa_map, id.vars = c("pointer"), measure.vars = c("OQ818324","OQ818318","PP711948",
                                                                                 "OR951334"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "OQ818324"] <- "R. madagascariensis teschovirus 2: OQ818324*"
long.sim_aa$accession[long.sim_aa$accession == "OQ818318"] <- "E. dupreanum teschovirus 1: OQ818318*"
long.sim_aa$accession[long.sim_aa$accession == "OR951334"] <- "Pteropodidae bat teschovirus: OR951334"
long.sim_aa$accession[long.sim_aa$accession == "PP711948"] <- "Rousettus bat picornavirus: PP711948"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("R. madagascariensis teschovirus 2: OQ818324*",
                                                                  "E. dupreanum teschovirus 1: OQ818318*",
                                                                  "Pteropodidae bat teschovirus: OR951334",
                                                                  "Rousettus bat picornavirus: PP711948"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## nucleotide
title<-expression(paste("Reference: R. madagascariensis teschovirus 1: OQ818323*"))

teschovirus_map_aa <- ggplot(long.sim_aa) + 
  annotate("rect", xmin=3368/3.055, xmax=3604/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=4961/3.055, xmax=5279/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Amino acid identity")+xlab("Genome position")+
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
  #scale_fill_distiller()+
  ggtitle(title)+
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

teschovirus_map_aa

#put gene map with PySimPlot
tescho_map_aa<-teschovirus_map_aa/map_tescho+plot_layout(nrow=2,  heights = c(1, 0.2))
tescho_map_aa

tescho_map_aa<-as.ggplot(tescho_map_aa)
tescho_map_aa




#Sapovirus
#Plot nucleotide
sapovirus_nt_map <- read.csv(file = "sapo_align_nt.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapovirus_nt_map)

#move to long
long.sim_nt <- melt(sapovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818319","PP712015",
                                                                               "PP712001","KX759623"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818319"] <- "E. dupreanum sapovirus 1: OQ818319*"
long.sim_nt$accession[long.sim_nt$accession == "PP712015"] <- "Rousettus bat calicivirus: PP712015"
long.sim_nt$accession[long.sim_nt$accession == "PP712001"] <- "Rousettus bat calicivirus: PP712001"
long.sim_nt$accession[long.sim_nt$accession == "KX759623"] <- "Bat sapovirus: KX759623"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum sapovirus 1: OQ818319*",
                                                                  "Rousettus bat calicivirus: PP712015",
                                                                  "Rousettus bat calicivirus: PP712001",
                                                                  "Bat sapovirus: KX759623"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Reference: E. dupreanum sapovirus 1: PP766459*"))


sapo_map_nt <- ggplot(long.sim_nt) + 
  annotate("rect", xmin=2000, xmax=2800, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5179, xmax=6800, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Nucleotide identity")+xlab("Genome position")+
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
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapo_map_nt

#put gene map with PySimPlot
sapo_nt<-sapo_map_nt/map_sapo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapo_nt

sapo_nt<-as.ggplot(sapo_nt)
sapo_nt

#amino acid
sapovirus_aa_map <- read.csv(file = "sapo_align_aa.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapovirus_aa_map)

#move to long
long.sim_aa <- melt(sapovirus_aa_map, id.vars = c("pointer"), measure.vars = c("OQ818319","PP712015",
                                                                               "PP712001","KX759623"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "OQ818319"] <- "E. dupreanum sapovirus 1: OQ818319*"
long.sim_aa$accession[long.sim_aa$accession == "PP712015"] <- "Rousettus bat calicivirus: PP712015"
long.sim_aa$accession[long.sim_aa$accession == "PP712001"] <- "Rousettus bat calicivirus: PP712001"
long.sim_aa$accession[long.sim_aa$accession == "KX759623"] <- "Bat sapovirus: KX759623"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("E. dupreanum sapovirus 1: OQ818319*",
                                                                  "Rousettus bat calicivirus: PP712015",
                                                                  "Rousettus bat calicivirus: PP712001",
                                                                  "Bat sapovirus: KX759623"))

#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

#plot nucleotide
title<-expression(paste("Reference: E. dupreanum sapovirus 1: PP766459*"))

sapovirus_map_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% Amino acid identity")+xlab("Genome position")+
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
  annotate("rect", xmin=2000/3.055, xmax=2800/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5179/3.055, xmax=6800/3.055, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  scale_color_manual(values=colzpalette) + 
  #scale_fill_distiller()+
  ggtitle(title)+
  coord_cartesian(ylim=c(0,1.02))+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_map_aa

#put gene map with PySimPlot
sapo_map_aa<-sapovirus_map_aa/map_sapo+plot_layout(nrow=2,  heights = c(1,0.2))
sapo_map_aa

sapo_map_aa<-as.ggplot(sapo_map_aa)
sapo_map_aa




#To plot figure 3
fig3<-plot_grid(hep_map_aa,kun_map_aa,bat_picorna_aa,
                 sapelo_map_aa,tescho_map_aa,sapo_map_aa,
                 ncol=3,
                 labels="AUTO", label_size = 23, align = "hv", axis="b")
fig3
#excluded kobuvirus from figure because it's so similar, and excluded mischivirus since it's so different


#To plot nucleotide supplementary
nucleotide_supp<-plot_grid(hepato_nt,kobu_nt,kunsagi_nt,mischi_nt,
                      bat_picorna_nt,sapelo_nt,tescho_nt,sapo_nt,
                      ncol=3,
                      labels="AUTO", label_size = 23, align = "hv", axis="b")
nucleotide_supp


#To plot amino acid supplementary for kobuvirus and mischivirus
amino_supp<-plot_grid(kobu_map_aa, mischi_map_aa,
                    ncol=2,
                    labels="AUTO", label_size = 23, align = "hv", axis="b")
amino_supp



# save figs

homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 

ggsave(file = paste0(homewd, "/final_figures/supplemental/Sfig4_nucleotide_pysimplot.pdf"),
       plot = nucleotide_supp,
       units="mm",  
       width=150, 
       height=110, 
       scale=4, 
       dpi=300)

ggsave(file = paste0(homewd, "/final_figures/Fig3_pysimplot.pdf"),
       plot = fig3,
       units="mm",  
       width=150, 
       height=80, 
       scale=4, 
       dpi=300)

ggsave(file = paste0(homewd, "/final_figures/supplemental/Sfig3_amino_pysimplot.pdf"),
       plot = amino_supp,
       units="mm",  
       width=100, 
       height=40, 
       scale=4, 
       dpi=300)



