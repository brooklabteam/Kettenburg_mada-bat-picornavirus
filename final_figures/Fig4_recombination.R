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
                                      "VP1","VP1 ", "2A", "2B", "2C", "3A", "3B",
                                      "3C", "3D", "NS1/NS2","Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Minor structural protein","3'UTR"))
#Load the peptide files
map_pep <- read.csv("bootscan_alignment_peptides.csv", header = T, stringsAsFactors = F)
map_pep$gene<-factor(map_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                              "VP1", "VP1 ","2A", "2B", "2C", "3A", "3B",
                                              "3C", "3D", "NS1/NS2","Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Minor structural protein","3'UTR"))
#Load the feature file in case its needed
map_feat <- read.csv("bootscan_alignment_features.csv", header = T, stringsAsFactors = F)
map_feat$gene<-factor(map_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                "VP1","VP1 ", "2A", "2B", "2C", "3A", "3B",
                                                "3C", "3D", "NS1/NS2","Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Minor structural protein","3'UTR"))

#Pick colors for genes
colz=c("5'UTR"="gold", "L"="royalblue","VP4"="paleturquoise3", "VP2"="skyblue1", "VP0"="royalblue4", "VP3"="steelblue1",
       "VP1"="cadetblue1", "2A"="orange1", "2B"="sienna2", "2C"="darkorange1", "3A"="palevioletred1", "3B"="plum",
       "3C"="rosybrown1", "3D"="pink3", "Helicase"="darkseagreen1","NS4"="darkolivegreen1","Vpg"="seagreen1","Pro-Pol"="palegreen1","NS1/NS2"="darkgreen",
       "Polyprotein"="azure3",
       "Minor structural protein"="black","3'UTR"="gold")

colz2=c("5'UTR"="grey", "L"="white","VP4"="white", "VP2"="white", "VP0"="white", "VP3"="white",
        "VP1"="white", "VP1 "="cornflowerblue", "2A"="cornflowerblue", "2B"="white", "2C"="white", "3A"="cornflowerblue", "3B"="white",
        "3C"="white", "3D"="white", "Helicase"="white","NS4"="cornflowerblue","Vpg"="white","Pro-Pol"="white","NS1/NS2"="white",
        "Polyprotein"="black",
        "Minor structural protein"="black","3'UTR"="grey")



#Plot map and BLAST together plots

#Subset map data by virus
map_hepatovirus<-subset(map,molecule=="Hepatovirus")
map_kobuvirus<-subset(map,molecule=="Kobuvirus")
map_sapovirus<-subset(map,molecule=="Sapovirus")
map_sapelovirus<-subset(map,molecule=="Sapelovirus")
map_teschovirus<-subset(map,molecule=="Teschovirus")

map_hepatovirus_pep<-subset(map_pep,molecule=="Hepatovirus")
map_sapovirus_pep<-subset(map_pep,molecule=="Sapovirus")
map_sapelovirus_pep<-subset(map_pep,molecule=="Sapelovirus")
map_teschovirus_pep<-subset(map_pep,molecule=="Teschovirus")

map_hepatovirus_feat<-subset(map_feat,molecule=="Hepatovirus")
map_sapovirus_feat<-subset(map_feat,molecule=="Sapovirus")
map_sapelovirus_feat<-subset(map_feat,molecule=="Sapelovirus")
map_teschovirus_feat<-subset(map_feat,molecule=="Teschovirus")

#gene maps
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
setwd("~/Desktop/developer/mada-bat-picornavirus/recombination/rdp_bootscan_csv_to_plot")

#colzpalette<-c("#F8766D","#C49A00","#53B400","#A58AFF","#00B6EB","darkorange1","#FB61D7")
colzpalette<-c("#3B9AB2","#EBCC2A","#F21A00")


#Hepatovirus
# Using NC_028366.1 as reference
hepato_map_nt <- read.csv(file = "hepato_rdp_NC028366view.csv", header = T, stringsAsFactors = F) #amino acid
head(hepato_map_nt)

#move to long
long.sim_nt <- melt(hepato_map_nt, id.vars = c("pointer"), measure.vars = c("NC_028981...PP766455","NC_028981...NC_028366",
                                                                            "PP766455...NC_028366"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_028981...PP766455"] <- "Unknown (inferred by Tupaia hepatovirus) (major parent) - E. dupreanum hepatovirus: PP766455* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "NC_028981...NC_028366"] <- "Unknown (inferred by Tupaia hepatovirus)  (major parent) - E. helvum hepatovirus M32Eidhel2010 (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "PP766455...NC_028366"] <- "E. dupreanum hepatovirus: PP766455* (minor parent) - E. helvum hepatovirus M32Eidhel2010 (recombinant)"


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

#5 tests sig from 1-790 then 7740 to 7917bp
#maybe sig also from 3495 to 3885bp
#NC_028981.1 used to infer unknown parent

#put gene map with bootscan
hepato_nt<-hepato_map_nt/map_hepato+plot_layout(nrow=2,  heights = c(1, 0.2))
hepato_nt

hepato_nt<-as.ggplot(hepato_nt)
hepato_nt

# Using PP766457 as reference
hepato_map_nt <- read.csv(file = "hepato_rdp_clade_PP766457view.csv", header = T, stringsAsFactors = F) #amino acid
head(hepato_map_nt)

#move to long
long.sim_nt <- melt(hepato_map_nt, id.vars = c("pointer"), measure.vars = c("NC_028366...KT452729","NC_028366...PP766457",
                                                                            "KT452729...PP766457"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_028366...KT452729"] <- "E.helvum hepatovirus M32Eidhel2010 (major parent) - R. landeri hepatovirus FO1AF48Rhilan2010 (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "NC_028366...PP766457"] <- "E.helvum hepatovirus M32Eidhel2010 (major parent) - E. dupreanum hepatovirus: PP766457* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "KT452729...PP766457"] <- "R. landeri hepatovirus FO1AF48Rhilan2010 (minor parent) - E. dupreanum hepatovirus: PP766457* (recombinant)"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E.helvum hepatovirus M32Eidhel2010 (major parent) - R. landeri hepatovirus FO1AF48Rhilan2010 (minor parent)",
                                                                  "E.helvum hepatovirus M32Eidhel2010 (major parent) - E. dupreanum hepatovirus: PP766457* (recombinant)",
                                                                  "R. landeri hepatovirus FO1AF48Rhilan2010 (minor parent) - E. dupreanum hepatovirus: PP766457* (recombinant)"))

long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Potential recombinant - E. dupreanum hepatovirus: PP766457*"))

hepato_map_clade_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  annotate("rect", xmin=1854, xmax=2077, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepato_map_clade_nt

#3 tests sig from 1854 to 2077nt, major is NC_028366 and minor is KT452729

#put gene map with bootscan
hepato_clade_nt<-hepato_map_clade_nt/map_hepato+plot_layout(nrow=2,  heights = c(1, 0.2))
hepato_clade_nt

hepato_clade_nt<-as.ggplot(hepato_clade_nt)
hepato_clade_nt




#Sapelovirus
#reference sequence OQ818321
sapelovirus_nt_map <- read.csv(file = "sapelo_rdp_clade_OQ818321view.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapelovirus_nt_map)

#move to long
long.sim_nt <- melt(sapelovirus_nt_map, id.vars = c("pointer"), measure.vars = c("PP711911...PP711943","PP711911...OQ818321",
                                                                                 "PP711943...OQ818321"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "PP711911...PP711943"] <- "Unknown (inferred by rousettus bat picornavirus 10A/Kenya/BAT22/2015) (major parent) - Eidolon bat picornavirus 6A/Kenya/BAT606/2015 (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "PP711911...OQ818321"] <- "Unknown (inferred by rousettus bat picornavirus 10A/Kenya/BAT22/2015) (major parent) - E. dupreanum sapelovirus 2: OQ818321* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "PP711943...OQ818321"] <- "Eidolon bat picornavirus 6A/Kenya/BAT606/2015 (minor parent) - E. dupreanum sapelovirus 2: OQ818321* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Unknown (inferred by rousettus bat picornavirus 10A/Kenya/BAT22/2015) (major parent) - Eidolon bat picornavirus 6A/Kenya/BAT606/2015 (minor parent)",
                                                                  "Unknown (inferred by rousettus bat picornavirus 10A/Kenya/BAT22/2015) (major parent) - E. dupreanum sapelovirus 2: OQ818321* (recombinant)",
                                                                  "Eidolon bat picornavirus 6A/Kenya/BAT606/2015 (minor parent) - E. dupreanum sapelovirus 2: OQ818321* (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Potential recombinant - E. dupreanum sapelovirus 2: OQ818321*"))

sapelo_map_clade_nt1 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  annotate("rect", xmin=336, xmax=957, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelo_map_clade_nt1

#5 tests sig, major parent unknown inferred by PP711911, minor parent PP711943 from 336-957nt

#put gene map with bootscan
sapelo_clade_nt1<-sapelo_map_clade_nt1/map_sapelo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapelo_clade_nt1

sapelo_clade_nt1<-as.ggplot(sapelo_clade_nt1)
sapelo_clade_nt1

#using NC_033820 as reference
sapelovirus_nt_map <- read.csv(file = "sapelo_rdp_clade_NC_033820view.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapelovirus_nt_map)

#move to long
long.sim_nt <- melt(sapelovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818320...PP711943","OQ818320...NC_033820",
                                                                                 "PP711943...NC_033820"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818320...PP711943"] <- "E. dupreanum sapelovirus 1: OQ818320* (major parent) - Eidolon bat picornavirus 6A/Kenya/BAT606/2015 (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818320...NC_033820"] <- "E. dupreanum sapelovirus 1: OQ818320* (major parent) - Bat sapelovirus Bat/CAM/Sap-P24/2013 (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "PP711943...NC_033820"] <- "Eidolon bat picornavirus 6A/Kenya/BAT606/2015 (minor parent) - Bat sapelovirus Bat/CAM/Sap-P24/2013 (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum sapelovirus 1: OQ818320* (major parent) - Eidolon bat picornavirus 6A/Kenya/BAT606/2015 (minor parent)",
                                                                  "E. dupreanum sapelovirus 1: OQ818320* (major parent) - Bat sapelovirus Bat/CAM/Sap-P24/2013 (recombinant)",
                                                                  "Eidolon bat picornavirus 6A/Kenya/BAT606/2015 (minor parent) - Bat sapelovirus Bat/CAM/Sap-P24/2013 (recombinant)"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## nucleotide
title<-expression(paste("Potential recombinant - Bat sapelovirus Bat/CAM/Sap-P24/2013"))

sapelo_map_clade_nt2 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  annotate("rect", xmin=3334, xmax=5410, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelo_map_clade_nt2

#5 tests sig, major parent OQ818320, minor parent PP711943 from 3334-5410nt

#put gene map with bootscan
sapelo_clade_nt2<-sapelo_map_clade_nt2/map_sapelo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapelo_clade_nt2

sapelo_clade_nt2<-as.ggplot(sapelo_clade_nt2)
sapelo_clade_nt2




#Sapovirus
#reference of E. helvum sapovirus
sapovirus_nt_map <- read.csv(file = "sapo_rdp_clade_ehelvumview.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapovirus_nt_map)

#move to long
long.sim_nt <- melt(sapovirus_nt_map, id.vars = c("pointer"), measure.vars = c("r_aegypticus_kenya...PP766459","r_aegypticus_kenya...e_helvum_cameroon",
                                                                               "PP766459...e_helvum_cameroon"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "r_aegypticus_kenya...PP766459"] <- "Unknown (inferred by r. aegypticus kenyan sapoviruses) (major parent) - E. dupreanum sapovirus 1: PP766459* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "r_aegypticus_kenya...e_helvum_cameroon"] <- "Unknown (inferred by r. aegypticus kenyan sapoviruses (major parent) - E. helvum sapovirus clade (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "PP766459...e_helvum_cameroon"] <- "E. dupreanum sapovirus 1: PP766459* (minor parent) - E. helvum sapovirus clade (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Unknown (inferred by r. aegypticus kenyan sapoviruses) (major parent) - E. dupreanum sapovirus 1: PP766459* (minor parent)",
                                                                  "Unknown (inferred by r. aegypticus kenyan sapoviruses (major parent) - E. helvum sapovirus clade (recombinant)",
                                                                  "E. dupreanum sapovirus 1: PP766459* (minor parent) - E. helvum sapovirus clade (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Potential recombinant - E. helvum sapovirus clade"))

sapo_map_clade_nt1 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  annotate("rect", xmin=3548, xmax=7519, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapo_map_clade_nt1

# major unknown inferred by r. aegypticus kenya, minor is PP766459, 4 tests sig from 3548 to 7519nt

#put gene map with bootscan
sapo_clade_nt1<-sapo_map_clade_nt1/map_sapo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapo_clade_nt1

sapo_clade_nt1<-as.ggplot(sapo_clade_nt1)
sapo_clade_nt1

#reference of OP963623
sapovirus_nt_map <- read.csv(file = "sapo_rdp_clade_OP963623view.csv", header = T, stringsAsFactors = F) #Amino acid
head(sapovirus_nt_map)

#move to long
long.sim_nt <- melt(sapovirus_nt_map, id.vars = c("pointer"), measure.vars = c("PP766459...OQ709197","PP766459...OP963623",
                                                                               "OQ709197...OP963623"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "PP766459...OQ709197"] <- "E. dupreanum sapovirus 1: PP766459* (major parent) - Bat sapovirus WD3 (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "PP766459...OP963623"] <- "E. dupreanum sapovirus 1: PP766459* (major parent) - Bat sapovirus BtSY2 (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ709197...OP963623"] <- "Bat sapovirus WD3 (minor parent) - Bat sapovirus BtSY2 (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("E. dupreanum sapovirus 1: PP766459* (major parent) - Bat sapovirus WD3 (minor parent)",
                                                                  "E. dupreanum sapovirus 1: PP766459* (major parent) - Bat sapovirus BtSY2 (recombinant)",
                                                                  "Bat sapovirus WD3 (minor parent) - Bat sapovirus BtSY2 (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Potential recombinant - Bat sapovirus BtSY2"))

sapo_map_clade_nt2 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  annotate("rect", xmin=7302, xmax=7504, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapo_map_clade_nt2

# major is PP766459, minor is OQ709197, 4 tests sig from 7302 to 7504

#put gene map with bootscan
sapo_clade_nt2<-sapo_map_clade_nt2/map_sapo+plot_layout(nrow=2,  heights = c(1, 0.2))
sapo_clade_nt2

sapo_clade_nt2<-as.ggplot(sapo_clade_nt2)
sapo_clade_nt2




#teschovirus
#reference of R. aegypticus africa
teschovirus_nt_map <- read.csv(file = "tescho_rdp_clade_raegypticusafricaview.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_nt_map)

#move to long
long.sim_nt <- melt(teschovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818323...OQ818324","OQ818323...r_aegypticus_africa",
                                                                                 "OQ818324...r_aegypticus_africa"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818323...OQ818324"] <- "R. madagascariensis teschovirus 1: OQ818323* (major parent) - Unknown (inferred by r. madagascariensis teschovirus 2: OQ818324*) (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818323...r_aegypticus_africa"] <- "R. madagascariensis teschovirus 1: OQ818323* (major parent) - R. aegypticus African teschovirus clade (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324...r_aegypticus_africa"] <- "Unknown (inferred by r. madagascariensis teschovirus 2: OQ818324*) (minor parent) - R. aegypticus African teschovirus clade (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. madagascariensis teschovirus 1: OQ818323* (major parent) - Unknown (inferred by r. madagascariensis teschovirus 2: OQ818324*) (minor parent)",
                                                                  "R. madagascariensis teschovirus 1: OQ818323* (major parent) - R. aegypticus African teschovirus clade (recombinant)",
                                                                  "Unknown (inferred by r. madagascariensis teschovirus 2: OQ818324*) (minor parent) - R. aegypticus African teschovirus clade (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Potential recombinant - R. aegypticus African teschovirus clade"))

tescho_map_clade_nt1 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  annotate("rect", xmin=7346, xmax=7558, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

tescho_map_clade_nt1

# major OQ818323 minor unknown inferred by OQ818324 6 tests sig from 7346 to 7558nt

#put gene map with bootscan
tescho_clade_nt1<-tescho_map_clade_nt1/map_tescho+plot_layout(nrow=2,  heights = c(1, 0.2))
tescho_clade_nt1

tescho_clade_nt1<-as.ggplot(tescho_clade_nt1)
tescho_clade_nt1

#reference of OQ818323
teschovirus_nt_map <- read.csv(file = "tescho_rdp_clade_OQ818323view.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_nt_map)

#move to long
long.sim_nt <- melt(teschovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818324...r_aegypticus_africa","OQ818324...OQ818323",
                                                                                 "r_aegypticus_africa...OQ818323"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818324...r_aegypticus_africa"] <- "R. madagascariensis teschovirus 2: OQ818324* (major parent) - R. aegypticus African teschovirus clade (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324...OQ818323"] <- "R. madagascariensis teschovirus 2: OQ818324* (major parent) - R. madagascariensis teschovirus 1: OQ818323* (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "r_aegypticus_africa...OQ818323"] <- "R. aegypticus African teschovirus clade (minor parent) - R. madagascariensis teschovirus 1: OQ818323* (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. madagascariensis teschovirus 2: OQ818324* (major parent) - R. aegypticus African teschovirus clade (minor parent)",
                                                                  "R. madagascariensis teschovirus 2: OQ818324* (major parent) - R. madagascariensis teschovirus 1: OQ818323* (recombinant)",
                                                                  "R. aegypticus African teschovirus clade (minor parent) - R. madagascariensis teschovirus 1: OQ818323* (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Potential recombinant - R. madagascariensis teschovirus 1: OQ818323*"))

tescho_map_clade_nt2 <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
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
  annotate("rect", xmin=63, xmax=1064, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

tescho_map_clade_nt2

# major parent OQ818324, minor OQ818323, 5 tests sig from 63-1064nt

#put gene map with bootscan
tescho_clade_nt2<-tescho_map_clade_nt2/map_tescho+plot_layout(nrow=2,  heights = c(1, 0.2))
tescho_clade_nt2

tescho_clade_nt2<-as.ggplot(tescho_clade_nt2)
tescho_clade_nt2

#reference of PP711934
teschovirus_nt_map <- read.csv(file = "tescho_rdp_PP711934view.csv", header = T, stringsAsFactors = F) #Amino acid
head(teschovirus_nt_map)

#move to long
long.sim_nt <- melt(teschovirus_nt_map, id.vars = c("pointer"), measure.vars = c("OQ818324...OQ818323","OQ818324...PP711934",
                                                                                 "OQ818323...PP711934"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818324...OQ818323"] <- "R. madagascariensis teschovirus 2: OQ818324* (major parent) - R. madagascariensis teschovirus 1: OQ818323* (minor parent)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324...PP711934"] <- "R. madagascariensis teschovirus 2: OQ818324* (major parent) - Rousettus bat picornavirus 29A/Kenya/BAT3/2015 (recombinant)"
long.sim_nt$accession[long.sim_nt$accession == "OQ818323...PP711934"] <- "R. madagascariensis teschovirus 1: OQ818323* (minor parent) - Rousettus bat picornavirus 29A/Kenya/BAT3/2015 (recombinant)"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("R. madagascariensis teschovirus 2: OQ818324* (major parent) - R. madagascariensis teschovirus 1: OQ818323* (minor parent)",
                                                                  "R. madagascariensis teschovirus 2: OQ818324* (major parent) - Rousettus bat picornavirus 29A/Kenya/BAT3/2015 (recombinant)",
                                                                  "R. madagascariensis teschovirus 1: OQ818323* (minor parent) - Rousettus bat picornavirus 29A/Kenya/BAT3/2015 (recombinant)"))

#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

#plot nucleotide
title<-expression(paste("Potential recombinant - Rousettus bat picornavirus 29A/Kenya/BAT3/2015"))

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
  geom_hline(yintercept=0.69, linetype='dashed', col = 'black')+
  annotate("rect", xmin=1, xmax=1104, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  annotate("rect", xmin=5395, xmax=7800, ymin=0, ymax=1, alpha=0.6,  fill="azure4")+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

tescho_map_nt1

#Major parent OQ818324, minor parent OQ818323, 5 tests sig from 1-1104 and 5395-7800nt

#put gene map with bootscan
tescho_nt1<-tescho_map_nt1/map_tescho+plot_layout(nrow=2,  heights = c(1, 0.2))
tescho_nt1

tescho_nt1<-as.ggplot(tescho_nt1)
tescho_nt1



#To plot figure 4
fig4<-plot_grid(hepato_clade_nt,sapelo_clade_nt2, tescho_nt1,
                ncol=3,
                labels="AUTO", label_size = 23, align = "hv", axis="b")
fig4


#To plot supplementary for everything else
recombination_supp<-plot_grid(hepato_nt,sapelo_clade_nt1, tescho_clade_nt1, tescho_clade_nt2, sapo_clade_nt1, sapo_clade_nt2,
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

ggsave(file = paste0(homewd, "/final_figures/supplemental/Sfig5_recombination.pdf"),
       plot = recombination_supp,
       units="mm",  
       width=150, 
       height=70, 
       scale=4, 
       dpi=300)


