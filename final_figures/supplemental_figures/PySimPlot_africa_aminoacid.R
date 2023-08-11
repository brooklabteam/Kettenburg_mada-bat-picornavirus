rm(list=ls())

library(ggplot2)
library(gggenes)t
library(cowplot)
library(gridExtra)
library(grid)
library(plyr)
library(dplyr)
library(reshape2)
library(cowplot)


#This is to make a figure of the representative PySimPlots with their matching genome plots


##Furst make the gene maps
setwd("~/Desktop/mada-bat-picornavirus/genome_annotation_and_characterization/genus_gene_maps")

#load files and make the genes into factor data
africa <- read.csv("africa_align_genes.csv", header = T, stringsAsFactors = F)
africa$gene<-factor(africa$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                            "VP1", "2A", "2B", "2C", "3A", "3B",
                                            "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                            "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                            "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#load the files and make the peptides into factor data
africa_pep <- read.csv("africa_align_peptides.csv", header = T, stringsAsFactors = F)
africa_pep$gene<-factor(africa_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                    "VP1", "2A", "2B", "2C", "3A", "3B",
                                                    "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                                    "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                    "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the files with gene and peptide markers for feature labels
africa_feat <- read.csv("africa_align_features.csv", header = T, stringsAsFactors = F)
africa_feat$gene<-factor(africa_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
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

#Subset africa data by virus
africa_felisavirus<-subset(africa, molecule == "Felisavirus")
africa_hepatovirus<-subset(africa, molecule == "Hepatovirus")
africa_kobuvirus<-subset(africa, molecule == "Kobuvirus")
africa_kunsagivirus<-subset(africa, molecule == "Kunsagivirus")
africa_mischivirus<-subset(africa, molecule == "Mischivirus")
africa_sapovirus<-subset(africa, molecule == "Sapovirus")
africa_sapelovirus<-subset(africa, molecule == "Sapelovirus")

africa_kobuvirus_pep<-subset(africa_pep, molecule == "Kobuvirus")
africa_mischivirus_pep<-subset(africa_pep, molecule == "Mischivirus")
africa_sapelovirus_pep<-subset(africa_pep, molecule == "Sapelovirus")

africa_felisavirus_feat<-subset(africa_feat, molecule == "Felisavirus")
africa_hepatovirus_feat<-subset(africa_feat, molecule == "Hepatovirus")
africa_kobuvirus_feat<-subset(africa_feat, molecule == "Kobuvirus")
africa_kunsagivirus_feat<-subset(africa_feat, molecule == "Kunsagivirus")
africa_mischivirus_feat<-subset(africa_feat, molecule == "Mischivirus")
africa_sapovirus_feat<-subset(africa_feat, molecule == "Sapovirus")
africa_sapelovirus_feat<-subset(africa_feat, molecule == "Sapelovirus")


#plot african gene maps
africa_felis<-ggplot(africa_felisavirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_text(data=africa_felisavirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,5805),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_felis


africa_hepato<-ggplot(africa_hepatovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_text(data=africa_hepatovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,6180),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_hepato


africa_kobu<-ggplot(africa_kobuvirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_subgene_arrow(data = africa_kobuvirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                            xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=africa_kobuvirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,8379),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_kobu



africa_kun<-ggplot(africa_kunsagivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_text(data=africa_kunsagivirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,7333),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_kun



africa_mischi<-ggplot(africa_mischivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_subgene_arrow(data = africa_mischivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=africa_mischivirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,8886),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_mischi


africa_sapo<-ggplot(africa_sapovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_text(data=africa_sapovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,7759),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_sapo


africa_sapelo<-ggplot(africa_sapelovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_subgene_arrow(data = africa_sapelovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                   xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=africa_sapelovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
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
africa_sapelo




##Now get the plots for the PySimPlot, just the ICTV_BLAST ones, there will be a saparate PDF file with the table of African bat picorna similarities
##Start with the African bat picornavirus comparisons
setwd("~/Desktop/mada-bat-picornavirus/PySimPlot/africa_bat_picorna_comparisons/africa_bat_pysimplot")


#Felisavirus
felisavirus_aa_africa <- read.csv(file = "felisavirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
head(felisavirus_aa_africa)

#move to long
long.sim_aa <- melt(felisavirus_aa_africa, id.vars = c("pointer"), measure.vars = c("KX644943","OQ818341"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "OQ818341"] <- "Eidolon dupreanum felisavirus OQ818341"
long.sim_aa$accession[long.sim_aa$accession == "KX644943"] <- "Eidolon helvum felisavirus KX644943"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon dupreanum felisavirus OQ818341", "Eidolon helvum felisavirus KX644943"
))

#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
title<-expression(paste("Reference: ",italic("Pteropus rufus felisavirus "), "OQ818335"))

felisavirus_africa_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
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
        plot.title = element_text(size = 12, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055),
                     labels = c(0,2000, 4000,6000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))
felisavirus_africa_aa

#put gene map with PySimPlot
felis_aa_africa<-felisavirus_africa_aa/africa_felis+plot_layout(nrow=2,  heights = c(2, 0.27))
felis_aa_africa

felis_aa_africa<-as.ggplot(felis_aa_africa)
felis_aa_africa



#Hepatovirus
hepatovirus_aa_africa <- read.csv(file = "hepatovirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
head(hepatovirus_aa_africa)

#move to long
long.sim_aa <- melt(hepatovirus_aa_africa, id.vars = c("pointer"), measure.vars = c("KT452712.1","KT452713.1","KT452729.1",
                                                                                    "KT452731.1","KT452732.1","KT452743.1",
                                                                                    "KT452744.1","NC_028366.1","NC_038313.1",
                                                                                    "NC_038316.1"))
unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "KT452712.1"] <- "Eidolon helvum hepatovirus KT452712"
long.sim_aa$accession[long.sim_aa$accession == "KT452713.1"] <- "Eidolon helvum hepatovirus KT452713"
long.sim_aa$accession[long.sim_aa$accession == "KT452729.1"] <- "Rhinolophus landeri hepatovirus KT452729"
long.sim_aa$accession[long.sim_aa$accession == "KT452731.1"] <- "Coleura afra hepatovirus KT452731"
long.sim_aa$accession[long.sim_aa$accession == "KT452732.1"] <- "Glauconycteris hepatovirus KT452732"
long.sim_aa$accession[long.sim_aa$accession == "KT452743.1"] <- "Miniopterus hepatovirus KT452743"
long.sim_aa$accession[long.sim_aa$accession == "KT452744.1"] <- "Miniopterus hepatovirus KT452744"
long.sim_aa$accession[long.sim_aa$accession == "NC_028366.1"] <- "Eidolon helvum hepatovirus NC_028366"
long.sim_aa$accession[long.sim_aa$accession == "NC_038313.1"] <- "Miniopterus hepatovirus NC_038313"
long.sim_aa$accession[long.sim_aa$accession == "NC_038316.1"] <- "Coleura afra hepatovirus NC_038316"


long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon helvum hepatovirus KT452712", 
                                                                  "Eidolon helvum hepatovirus KT452713",
                                                                  "Rhinolophus landeri hepatovirus KT452729",
                                                                  "Coleura afra hepatovirus KT452731",
                                                                  "Glauconycteris hepatovirus KT452732",
                                                                  "Miniopterus hepatovirus KT452743",
                                                                  "Miniopterus hepatovirus KT452744",
                                                                  "Eidolon helvum hepatovirus NC_028366",
                                                                  "Miniopterus hepatovirus NC_038313",
                                                                  "Coleura afra hepatovirus NC_038316"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
title<-expression(paste("Reference: ",italic("Eidolon dupreanum hepatovirus "), "OQ818337"))

hepatovirus_africa_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 12, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055),
                     labels = c(0,2000, 4000,6000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_africa_aa

#put gene map with PySimPlot
hep_africa_aa<-hepatovirus_africa_aa/africa_hepato+plot_layout(nrow=2,  heights = c(2, 0.27))
hep_africa_aa

hep_africa_aa<-as.ggplot(hep_africa_aa)
hep_africa_aa




#Kobuvirus
kobuvirus_aa_africa <- read.csv(file = "kobuvirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
head(kobuvirus_aa_africa)

#move to long
long.sim_aa <- melt(kobuvirus_aa_africa, id.vars = c("pointer"), measure.vars = c("OP287812.1","JX885611.1"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "OP287812.1"] <- "Eidolon dupreanum kobuvirus OP287812"
long.sim_aa$accession[long.sim_aa$accession == "JX885611.1"] <- "Eidolon helvum kobuvirus JX885611"


long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon dupreanum kobuvirus OP287812", 
                                                                  "Eidolon helvum kobuvirus JX885611"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
title<-expression(paste("Reference: ",italic("Eidolon dupreanum kobuvirus "), "OQ818322"))

kobuvirus_africa_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
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
        plot.title = element_text(size = 12, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobuvirus_africa_aa

#put gene map with PySimPlot
kobu_africa_aa<-kobuvirus_africa_aa/africa_kobu+plot_layout(nrow=2,  heights = c(2, 0.27))
kobu_africa_aa

kobu_africa_aa<-as.ggplot(kobu_africa_aa)
kobu_africa_aa



#Kunsagivirus
kunsagivirus_aa_africa <- read.csv(file = "kunsagivirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
head(kunsagivirus_aa_africa)

#move to long
long.sim_aa <- melt(kunsagivirus_aa_africa, id.vars = c("pointer"), measure.vars = c("NC_033818.1"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "NC_033818.1"] <- "Eidolon dupreanum helvum NC_033818"


long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon helvum kunsagivirus NC_033818"))

#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

colz2=c("Eidolon helvum kunsagivirus NC_033818"="black")

## amino acid
title<-expression(paste("Reference: ",italic("Eidolon dupreanum kunsagivirus "), "OQ818317"))

kunsagivirus_africa_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
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
        plot.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(values=colz2) + 
  #scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagivirus_africa_aa

#put gene map with PySimPlot
kun_africa_aa<-kunsagivirus_africa_aa/africa_kun+plot_layout(nrow=2,  heights = c(2, 0.27))
kun_africa_aa

kun_africa_aa<-as.ggplot(kun_africa_aa)
kun_africa_aa



#mischivirus
mischivirus_aa_africa <- read.csv(file = "mischivirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
head(mischivirus_aa_africa)

#move to long
long.sim_aa <- melt(mischivirus_aa_africa, id.vars = c("pointer"), measure.vars = c("MG888045.1","NC_026470.1"))

unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "MG888045.1"] <- "Miniopterus mischivirus MG888045"
long.sim_aa$accession[long.sim_aa$accession == "NC_026470.1"] <- "Macronycteris gigas mischivirus NC_026470"


long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Miniopterus mischivirus MG888045", 
                                                                  "Macronycteris gigas mischivirus NC_026470"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
title<-expression(paste("Reference: ",italic("Pteropus rufus mischivirus "), "OQ818316"))

mischivirus_africa_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
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
        plot.title = element_text(size = 12, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_africa_aa

#put gene map with PySimPlot
mischi_africa_aa<-mischivirus_africa_aa/africa_mischi+plot_layout(nrow=2,  heights = c(2, 0.27))
mischi_africa_aa

mischi_africa_aa<-as.ggplot(mischi_africa_aa)
mischi_africa_aa



#Sapelovirus
sapelovirus_aa_africa <- read.csv(file = "sapelovirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
head(sapelovirus_aa_africa)

#move to long
long.sim_aa <- melt(sapelovirus_aa_africa, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329",
                                                                                    "OQ818342","OQ818343",
                                                                                    "OQ818344","NC_033820.1"))
unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_aa$accession[long.sim_aa$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_aa$accession[long.sim_aa$accession == "OQ818342"] <- "Eidolon dupreanum sapelovirus OQ818342"
long.sim_aa$accession[long.sim_aa$accession == "OQ818343"] <- "Eidolon dupreanum sapelovirus OQ818343"
long.sim_aa$accession[long.sim_aa$accession == "OQ818344"] <- "Eidolon dupreanum sapelovirus OQ818344"
long.sim_aa$accession[long.sim_aa$accession == "NC_033820.1"] <- "Eidolon helvum sapelovirus NC_033820"


long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321", 
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Eidolon dupreanum sapelovirus OQ818342",
                                                                  "Eidolon dupreanum sapelovirus OQ818343",
                                                                  "Eidolon dupreanum sapelovirus OQ818344",
                                                                  "Eidolon helvum sapelovirus NC_033820"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapelovirus "), "OQ818320"))

sapelovirus_africa_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
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
        plot.title = element_text(size = 12, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_africa_aa

#put gene map with PySimPlot
sapelo_africa_aa<-sapelovirus_africa_aa/africa_sapelo+plot_layout(nrow=2,  heights = c(2, 0.27))
sapelo_africa_aa

sapelo_africa_aa<-as.ggplot(sapelo_africa_aa)
sapelo_africa_aa





#Sapovirus
sapovirus_aa_africa <- read.csv(file = "sapovirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
head(sapovirus_aa_africa)

#move to long
long.sim_aa <- melt(sapovirus_aa_africa, id.vars = c("pointer"), measure.vars =c("OQ818340","OQ818345","OQ818347","OQ818348",
                                                                                 "KX759618.1","KX759619.1","KX759621.1",
                                                                                 "KX759622.1","KX759623.1","KX759624.1",
                                                                                 "NC_033776.1"))
unique(long.sim_aa$variable)

long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_aa$accession[long.sim_aa$accession == "OQ818340"] <- "Eidolon dupreanum sapovirus OQ818340"
long.sim_aa$accession[long.sim_aa$accession == "OQ818345"] <- "Rousettus madagascariensis sapovirus OQ818345"
long.sim_aa$accession[long.sim_aa$accession == "OQ818347"] <- "Rousettus madagascariensis sapovirus OQ818347"
long.sim_aa$accession[long.sim_aa$accession == "OQ818348"] <- "Rousettus madagascariensis sapovirus OQ818348"
long.sim_aa$accession[long.sim_aa$accession == "KX759618.1"] <- "Eidolon helvum sapovirus KX759618"
long.sim_aa$accession[long.sim_aa$accession == "KX759619.1"] <- "Eidolon helvum sapovirus KX759619"
long.sim_aa$accession[long.sim_aa$accession == "KX759621.1"] <- "Eidolon helvum sapovirus KX759621"
long.sim_aa$accession[long.sim_aa$accession == "KX759622.1"] <- "Eidolon helvum sapovirus KX759622"
long.sim_aa$accession[long.sim_aa$accession == "KX759623.1"] <- "Eidolon helvum sapovirus KX759623"
long.sim_aa$accession[long.sim_aa$accession == "KX759624.1"] <- "Eidolon helvum sapovirus KX759624"
long.sim_aa$accession[long.sim_aa$accession == "NC_033776.1"] <- "Eidolon helvum sapovirus NC_033776"

long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon dupreanum sapovirus OQ818340",
                                                                  "Rousettus madagascariensis sapovirus OQ818345", 
                                                                  "Rousettus madagascariensis sapovirus OQ818347",
                                                                  "Rousettus madagascariensis sapovirus OQ818348",
                                                                  "Eidolon helvum sapovirus KX759618",
                                                                  "Eidolon helvum sapovirus KX759619",
                                                                  "Eidolon helvum sapovirus KX759621",
                                                                  "Eidolon helvum sapovirus KX759622",
                                                                  "Eidolon helvum sapovirus KX759623",
                                                                  "Eidolon helvum sapovirus NC_038313",
                                                                  "Eidolon helvum sapovirus NC_033776"))
#and plot
long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapovirus "), "OQ818319"))

sapovirus_africa_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+xlab("Genome position")+
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
        plot.title = element_text(size = 12, face = "bold")) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055),
                     labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  #scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_africa_aa

#put gene map with PySimPlot
sapo_africa_aa<-sapovirus_africa_aa/africa_sapo+plot_layout(nrow=2,  heights = c(2, 0.27))
sapo_africa_aa

sapo_africa_aa<-as.ggplot(sapo_africa_aa)
sapo_africa_aa



##Now put the whole figure together
supp_africa_pysimplot<-plot_grid(felis_aa_africa,hep_africa_aa,kobu_africa_aa,kun_africa_aa,
                mischi_africa_aa,sapelo_africa_aa,sapo_africa_aa,
                ncol=3,
                labels=c("A","","","","","",""))
supp_africa_pysimplot

