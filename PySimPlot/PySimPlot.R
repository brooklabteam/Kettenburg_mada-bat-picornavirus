rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)

homewd="/Users/gkettenburg/Desktop/mada-bat-picornavirus/PySimPlot"

colz=c("5'UTR"="gold", "L"="royalblue","VP4"="paleturquoise3", "VP2"="skyblue1", "VP0"="royalblue4", "VP3"="steelblue1",
       "VP1"="cadetblue1", "2A"="palevioletred1", "2B"="red4", "2C"="palevioletred3", "3A"="tomato2", "3B"="plum",
       "3C"="rosybrown1", "3D"="pink2", 
       "Polyprotein"="azure3","Putative polyprotein"="mediumorchid1", "Non-structural polyprotein"="mediumorchid4", 
       "Putative minor structural protein" ="slateblue3", "Structural polyprotein"="lightgoldenrod1",
       "Hypothetical protein"="darkslategrey", "Similar to structural polyprotein"="mediumpurple1", "Similar to putative polyprotein"="mediumpurple3", 
       "Similar to polyprotein"="mediumpurple4","3'UTR"="yellow")


##Start with the African bat picornavirus comparisons
setwd("~/Desktop/mada-bat-picornavirus/PySimPlot/africa_bat_picorna_comparisons/africa_bat_pysimplot")


#Felisavirus
felisavirus_aa_africa <- read.csv(file = "felisavirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
felisavirus_nt_africa <- read.csv(file = "felisavirus_nt_africa.csv", header = T, stringsAsFactors = F) #animo acid
head(felisavirus_nt_africa)
head(felisavirus_aa_africa)

#move to long
long.sim_nt <- melt(felisavirus_nt_africa, id.vars = c("pointer"), measure.vars = c("OQ818341","KX644943"))
long.sim_aa <- melt(felisavirus_aa_africa, id.vars = c("pointer"), measure.vars = c("KX644943","OQ818341"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818341"] <- "Eidolon dupreanum felisavirus OQ818341"
long.sim_nt$accession[long.sim_nt$accession == "KX644943"] <- "Eidolon helvum felisavirus KX644943"
long.sim_aa$accession[long.sim_aa$accession == "OQ818341"] <- "Eidolon dupreanum felisavirus OQ818341"
long.sim_aa$accession[long.sim_aa$accession == "KX644943"] <- "Eidolon helvum felisavirus KX644943"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum felisavirus OQ818341", "Eidolon helvum felisavirus KX644943" 
                                                            ))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon dupreanum felisavirus OQ818341", "Eidolon helvum felisavirus KX644943"
                                                            ))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

#Keep the script to make a banner in case I need it again
# genome_nt <- data.frame(position = c(1,1251,
#                                         1252,1537,
#                                         1538,1761,
#                                         1762,2566,
#                                         2567,3270,
#                                         3274,4198,
#                                         4199,4294,
#                                         4295,4817,
#                                         4818,5979,
#                                         5890,6241,
#                                         6262,6324,
#                                         6325,6961,
#                                         6962,8392,
#                                         8396, 8572), 
#                            gene = rep(c("5'UTR", "L peptide", "VP4 peptide", "VP2 peptide", "VP3 peptide","VP1 peptide", "2A peptide",
#                                         "2B peptide", "2C peptide", "3A peptide", "3B peptide", "3C peptide", "3D peptide", "3'UTR"), each=2))
# 
# genome_poly <- data.frame(position = c(1,1251,
#                                      1270, 8395,
#                                      8396, 8572), 
#                         gene = rep(c("5'UTR", "Polyprotein", "3'UTR"), each=2))
# 
# genome_nt$gene <- factor(genome_nt$gene, levels = c("5'UTR", "L peptide","VP4 peptide", "VP2 peptide", "VP0 peptide", "VP3 peptide",
#                                                     "VP1 peptide", "2A peptide", "2B peptide", "2C peptide", "3A peptide", "3B peptide",
#                                                     "3C peptide", "3D peptide", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein", 
#                                                     "Hypothetical protein", "Similar to structural polyprotein", 
#                                                     "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))

colz2= c("Eidolon dupreanum felisavirus OQ818341"="black", "Eidolon helvum felisavirus KX644943"="grey")

## amino acid
felisavirus_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Pteropus rufus felisavirus OQ818335")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

felisavirus_aa

## nucleotide
felisavirus_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

felisavirus_nt

library(cowplot)
felisavirus<-plot_grid(
  felisavirus_aa, felisavirus_nt,
  labels = NULL, ncol = 1
)
felisavirus



#Hepatovirus
hepatovirus_aa_africa <- read.csv(file = "hepatovirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
hepatovirus_nt_africa <- read.csv(file = "hepatovirus_nt_africa.csv", header = T, stringsAsFactors = F) #animo acid
head(hepatovirus_nt_africa)
head(hepatovirus_aa_africa)

#move to long
long.sim_nt <- melt(hepatovirus_nt_africa, id.vars = c("pointer"), measure.vars = c("KT452712.1","KT452713.1","KT452729.1",
                                                                                    "KT452731.1","KT452732.1","KT452743.1",
                                                                                    "KT452744.1","NC_028366.1","NC_038313.1",
                                                                                    "NC_038316.1"))
long.sim_aa <- melt(hepatovirus_aa_africa, id.vars = c("pointer"), measure.vars = c("KT452712.1","KT452713.1","KT452729.1",
                                                                                    "KT452731.1","KT452732.1","KT452743.1",
                                                                                    "KT452744.1","NC_028366.1","NC_038313.1",
                                                                                    "NC_038316.1"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "KT452712.1"] <- "Eidolon helvum hepatovirus KT452712"
long.sim_nt$accession[long.sim_nt$accession == "KT452713.1"] <- "Eidolon helvum hepatovirus KT452713"
long.sim_nt$accession[long.sim_nt$accession == "KT452729.1"] <- "Rhinolophus landeri hepatovirus KT452729"
long.sim_nt$accession[long.sim_nt$accession == "KT452731.1"] <- "Coleura afra hepatovirus KT452731"
long.sim_nt$accession[long.sim_nt$accession == "KT452732.1"] <- "Glauconycteris hepatovirus KT452732"
long.sim_nt$accession[long.sim_nt$accession == "KT452743.1"] <- "Miniopterus hepatovirus KT452743"
long.sim_nt$accession[long.sim_nt$accession == "KT452744.1"] <- "Miniopterus hepatovirus KT452744"
long.sim_nt$accession[long.sim_nt$accession == "NC_028366.1"] <- "Eidolon helvum hepatovirus NC_028366"
long.sim_nt$accession[long.sim_nt$accession == "NC_038313.1"] <- "Miniopterus hepatovirus NC_038313"
long.sim_nt$accession[long.sim_nt$accession == "NC_038316.1"] <- "Coleura afra hepatovirus NC_038316"

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


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon helvum hepatovirus KT452712", 
                                                                  "Eidolon helvum hepatovirus KT452713",
                                                                  "Rhinolophus landeri hepatovirus KT452729",
                                                                  "Coleura afra hepatovirus KT452731",
                                                                  "Glauconycteris hepatovirus KT452732",
                                                                  "Miniopterus hepatovirus KT452743",
                                                                  "Miniopterus hepatovirus KT452744",
                                                                  "Eidolon helvum hepatovirus NC_028366",
                                                                  "Miniopterus hepatovirus NC_038313",
                                                                  "Coleura afra hepatovirus NC_038316"))
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
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
hepatovirus_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_fill_distiller()+
  #scale_color_manual(values=colz2) + 
  ggtitle("Reference: Eidolon dupreanum hepatovirus OQ818337")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_aa

## nucleotide
hepatovirus_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_nt

library(cowplot)
hepatovirus<-plot_grid(
  hepatovirus_aa, hepatovirus_nt,
  labels = NULL, ncol = 1
)
hepatovirus



#Kobuvirus
kobuvirus_aa_africa <- read.csv(file = "kobuvirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
kobuvirus_nt_africa <- read.csv(file = "kobuvirus_nt_africa.csv", header = T, stringsAsFactors = F) #animo acid
head(kobuvirus_nt_africa)
head(kobuvirus_aa_africa)

#move to long
long.sim_nt <- melt(kobuvirus_nt_africa, id.vars = c("pointer"), measure.vars = c("OP287812.1","JX885611.1"))
long.sim_aa <- melt(kobuvirus_aa_africa, id.vars = c("pointer"), measure.vars = c("OP287812.1","JX885611.1"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OP287812.1"] <- "Eidolon dupreanum kobuvirus OP287812"
long.sim_nt$accession[long.sim_nt$accession == "JX885611.1"] <- "Eidolon helvum kobuvirus JX885611"

long.sim_aa$accession[long.sim_aa$accession == "OP287812.1"] <- "Eidolon dupreanum kobuvirus OP287812"
long.sim_aa$accession[long.sim_aa$accession == "JX885611.1"] <- "Eidolon helvum kobuvirus JX885611"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum kobuvirus OP287812", 
                                                                  "Eidolon helvum kobuvirus JX885611"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon dupreanum kobuvirus OP287812", 
                                                                  "Eidolon helvum kobuvirus JX885611"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
kobuvirus_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_fill_distiller()+
  #scale_color_manual(values=colz2) + 
  ggtitle("Reference: Eidolon dupreanum kobuvirus OQ818322")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobuvirus_aa

## nucleotide
kobuvirus_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobuvirus_nt

library(cowplot)
kobuvirus<-plot_grid(
  kobuvirus_aa, kobuvirus_nt,
  labels = NULL, ncol = 1
)
kobuvirus



#Kunsagivirus
kunsagivirus_aa_africa <- read.csv(file = "kunsagivirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
kunsagivirus_nt_africa <- read.csv(file = "kunsagivirus_nt_africa.csv", header = T, stringsAsFactors = F) #animo acid
head(kunsagivirus_nt_africa)
head(kunsagivirus_aa_africa)

#move to long
long.sim_nt <- melt(kunsagivirus_nt_africa, id.vars = c("pointer"), measure.vars = c("NC_033818.1"))
long.sim_aa <- melt(kunsagivirus_aa_africa, id.vars = c("pointer"), measure.vars = c("NC_033818.1"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_033818.1"] <- "Eidolon dupreanum helvum NC_033818"

long.sim_aa$accession[long.sim_aa$accession == "NC_033818.1"] <- "Eidolon dupreanum helvum NC_033818"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon helvum kunsagivirus NC_033818")) 
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon helvum kunsagivirus NC_033818"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

colz2=c("Eidolon helvum kunsagivirus NC_033818"="black")

## amino acid
kunsagivirus_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_fill_distiller()+
  scale_color_manual(values=colz2) + 
  ggtitle("Reference: Eidolon dupreanum kunsagivirus OQ818317")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagivirus_aa

## nucleotide
kunsagivirus_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_color_manual(values=colz2) + 
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(-0.1,1), expand=c(0,0))

kunsagivirus_nt

library(cowplot)
kunsagivirus<-plot_grid(
  kunsagivirus_aa, kunsagivirus_nt,
  labels = NULL, ncol = 1
)
kunsagivirus


#mischivirus
mischivirus_aa_africa <- read.csv(file = "mischivirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
mischivirus_nt_africa <- read.csv(file = "mischivirus_nt_africa.csv", header = T, stringsAsFactors = F) #animo acid
head(mischivirus_nt_africa)
head(mischivirus_aa_africa)

#move to long
long.sim_nt <- melt(mischivirus_nt_africa, id.vars = c("pointer"), measure.vars = c("MG888045.1","NC_026470.1"))
long.sim_aa <- melt(mischivirus_aa_africa, id.vars = c("pointer"), measure.vars = c("MG888045.1","NC_026470.1"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "MG888045.1"] <- "Miniopterus mischivirus MG888045"
long.sim_nt$accession[long.sim_nt$accession == "NC_026470.1"] <- "Macronycteris gigas mischivirus NC_026470"

long.sim_aa$accession[long.sim_aa$accession == "MG888045.1"] <- "Miniopterus mischivirus MG888045"
long.sim_aa$accession[long.sim_aa$accession == "NC_026470.1"] <- "Macronycteris gigas mischivirus NC_026470"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Miniopterus mischivirus MG888045", 
                                                                  "Macronycteris gigas mischivirus NC_026470"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Miniopterus mischivirus MG888045", 
                                                                  "Macronycteris gigas mischivirus NC_026470"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
mischivirus_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_fill_distiller()+
  #scale_color_manual(values=colz2) + 
  ggtitle("Reference: Pteropus rufus mischivirus OQ818316")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_aa

## nucleotide
mischivirus_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_nt

library(cowplot)
mischivirus<-plot_grid(
  mischivirus_aa, mischivirus_nt,
  labels = NULL, ncol = 1
)
mischivirus


#Sapelovirus
sapelovirus_aa_africa <- read.csv(file = "sapelovirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
sapelovirus_nt_africa <- read.csv(file = "sapelovirus_nt_africa.csv", header = T, stringsAsFactors = F) #animo acid
head(sapelovirus_nt_africa)
head(sapelovirus_aa_africa)

#move to long
long.sim_nt <- melt(sapelovirus_nt_africa, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329",
                                                                                    "OQ818342","OQ818343",
                                                                                    "OQ818344","NC_033820.1"))
long.sim_aa <- melt(sapelovirus_aa_africa, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329",
                                                                                    "OQ818342","OQ818343",
                                                                                    "OQ818344","NC_033820.1"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_nt$accession[long.sim_nt$accession == "OQ818342"] <- "Eidolon dupreanum sapelovirus OQ818342"
long.sim_nt$accession[long.sim_nt$accession == "OQ818343"] <- "Eidolon dupreanum sapelovirus OQ818343"
long.sim_nt$accession[long.sim_nt$accession == "OQ818344"] <- "Eidolon dupreanum sapelovirus OQ818344"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820.1"] <- "Eidolon helvum sapelovirus NC_033820"

long.sim_aa$accession[long.sim_aa$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_aa$accession[long.sim_aa$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_aa$accession[long.sim_aa$accession == "OQ818342"] <- "Eidolon dupreanum sapelovirus OQ818342"
long.sim_aa$accession[long.sim_aa$accession == "OQ818343"] <- "Eidolon dupreanum sapelovirus OQ818343"
long.sim_aa$accession[long.sim_aa$accession == "OQ818344"] <- "Eidolon dupreanum sapelovirus OQ818344"
long.sim_aa$accession[long.sim_aa$accession == "NC_033820.1"] <- "Eidolon helvum sapelovirus NC_033820"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321", 
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Eidolon dupreanum sapelovirus OQ818342",
                                                                  "Eidolon dupreanum sapelovirus OQ818343",
                                                                  "Eidolon dupreanum sapelovirus OQ818344",
                                                                  "Eidolon helvum sapelovirus NC_033820"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321", 
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Eidolon dupreanum sapelovirus OQ818342",
                                                                  "Eidolon dupreanum sapelovirus OQ818343",
                                                                  "Eidolon dupreanum sapelovirus OQ818344",
                                                                  "Eidolon helvum sapelovirus NC_033820"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
sapelovirus_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_fill_distiller()+
  #scale_color_manual(values=colz2) + 
  ggtitle("Reference: Eidolon dupreanum sapelovirus OQ818320")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_aa

## nucleotide
sapelovirus_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_nt

library(cowplot)
sapelovirus<-plot_grid(
  sapelovirus_aa, sapelovirus_nt,
  labels = NULL, ncol = 1
)
sapelovirus


#Sapovirus
sapovirus_aa_africa <- read.csv(file = "sapovirus_aa_africa.csv", header = T, stringsAsFactors = F) #nucleotide
sapovirus_nt_africa <- read.csv(file = "sapovirus_nt_africa.csv", header = T, stringsAsFactors = F) #animo acid
head(sapovirus_nt_africa)
head(sapovirus_aa_africa)

#move to long
long.sim_nt <- melt(sapovirus_nt_africa, id.vars = c("pointer"), measure.vars = c("OQ818340","OQ818345","OQ818347","OQ818348",
                                                                                    "KX759618.1","KX759619.1","KX759621.1",
                                                                                    "KX759622.1","KX759623.1","KX759624.1",
                                                                                    "NC_033776.1"))
long.sim_aa <- melt(sapovirus_aa_africa, id.vars = c("pointer"), measure.vars =c("OQ818340","OQ818345","OQ818347","OQ818348",
                                                                                 "KX759618.1","KX759619.1","KX759621.1",
                                                                                 "KX759622.1","KX759623.1","KX759624.1",
                                                                                 "NC_033776.1"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818340"] <- "Eidolon dupreanum sapovirus OQ818340"
long.sim_nt$accession[long.sim_nt$accession == "OQ818345"] <- "Rousettus madagascariensis sapovirus OQ818345"
long.sim_nt$accession[long.sim_nt$accession == "OQ818347"] <- "Rousettus madagascariensis sapovirus OQ818347"
long.sim_nt$accession[long.sim_nt$accession == "OQ818348"] <- "Rousettus madagascariensis sapovirus OQ818348"
long.sim_nt$accession[long.sim_nt$accession == "KX759618.1"] <- "Eidolon helvum sapovirus KX759618"
long.sim_nt$accession[long.sim_nt$accession == "KX759619.1"] <- "Eidolon helvum sapovirus KX759619"
long.sim_nt$accession[long.sim_nt$accession == "KX759621.1"] <- "Eidolon helvum sapovirus KX759621"
long.sim_nt$accession[long.sim_nt$accession == "KX759622.1"] <- "Eidolon helvum sapovirus KX759622"
long.sim_nt$accession[long.sim_nt$accession == "KX759623.1"] <- "Eidolon helvum sapovirus KX759623"
long.sim_nt$accession[long.sim_nt$accession == "KX759624.1"] <- "Eidolon helvum sapovirus KX759624"
long.sim_nt$accession[long.sim_nt$accession == "NC_033776.1"] <- "Eidolon helvum sapovirus NC_033776"

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


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum sapovirus OQ818340",
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
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
sapovirus_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_fill_distiller()+
  #scale_color_manual(values=colz2) + 
  ggtitle("Reference: Eidolon dupreanum sapovirus OQ818319")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_aa

## nucleotide
sapovirus_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_nt

library(cowplot)
sapovirus<-plot_grid(
  sapovirus_aa, sapovirus_nt,
  labels = NULL, ncol = 1
)
sapovirus





##Now do ICTV+BLAST comparisons
setwd("~/Desktop/mada-bat-picornavirus/PySimPlot/ICTV_BLAST_pysimplot")


#Bat picornavirus 
bat_picorna_ictv_nt_full <- read.csv(file = "bat_picorna_ictv_nt_full_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
bat_picorna_ictv_aa_full <- read.csv(file = "bat_picorna_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(bat_picorna_ictv_nt_full)
head(bat_picorna_ictv_aa_full)

#move to long
long.sim_nt <- melt(bat_picorna_ictv_nt_full, id.vars = c("pointer"), measure.vars = c("OQ818328","HQ595340","HQ595342","HQ595344"))
long.sim_aa <- melt(bat_picorna_ictv_aa_full, id.vars = c("pointer"), measure.vars = c("OQ818328","HQ595340","HQ595342","HQ595344"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818328"] <- "Rousettus madagascariensis picornavirus OQ818328"
long.sim_nt$accession[long.sim_nt$accession == "HQ595340"] <- "Bat picornavirus 1 HQ595340"
long.sim_nt$accession[long.sim_nt$accession == "HQ595342"] <- "Bat picornavirus 2 HQ595342"
long.sim_nt$accession[long.sim_nt$accession == "HQ595344"] <- "Bat picornavirus 3 HQ595344"

long.sim_aa$accession[long.sim_aa$accession == "OQ818328"] <- "Rousettus madagascariensis picornavirus OQ818328"
long.sim_aa$accession[long.sim_aa$accession == "HQ595340"] <- "Bat picornavirus 1 HQ595340"
long.sim_aa$accession[long.sim_aa$accession == "HQ595342"] <- "Bat picornavirus 2 HQ595342"
long.sim_aa$accession[long.sim_aa$accession == "HQ595344"] <- "Bat picornavirus 3 HQ595344"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rousettus madagascariensis picornavirus OQ818328", "Bat picornavirus 1 HQ595340",
                                                                  "Bat picornavirus 2 HQ595342","Bat picornavirus 3 HQ595344"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Rousettus madagascariensis picornavirus OQ818328", "Bat picornavirus 1 HQ595340",
                                                                  "Bat picornavirus 2 HQ595342","Bat picornavirus 3 HQ595344"))
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
batpicorna_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Rousettus madagascariensis picornavirus OQ818325")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_ictv_aa

## nucleotide
batpicorna_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_ictv_nt

library(cowplot)
bat_picorna_ictv<-plot_grid(
  batpicorna_ictv_aa, batpicorna_ictv_nt,
  labels = NULL, ncol = 1
)
bat_picorna_ictv



#Cheravirus RNA2
cheravirus_ictv_nt_full <- read.csv(file = "chera_ictv_nt_partial_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
cheravirus_ictv_aa_full <- read.csv(file = "chera_ictv_aa_partial_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(cheravirus_ictv_nt_full)
head(cheravirus_ictv_aa_full)

#move to long
long.sim_nt <- melt(cheravirus_ictv_nt_full, id.vars = c("pointer"), measure.vars = c("AB030941","AJ621358","KT692953","DQ143875","MK153132","MW582786"))
long.sim_aa <- melt(cheravirus_ictv_aa_full, id.vars = c("pointer"), measure.vars = c("AB030941","AJ621358","KT692953","DQ143875","MK153132","MW582786"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "AB030941"] <- "Cheravirus mali AB030941"
long.sim_nt$accession[long.sim_nt$accession == "AJ621358"] <- "Cheravirus avii AJ621358"
long.sim_nt$accession[long.sim_nt$accession == "KT692953"] <- "Cheravirus ribis KT692953"
long.sim_nt$accession[long.sim_nt$accession == "DQ143875"] <- "Cheravirus pruni DQ143875"
long.sim_nt$accession[long.sim_nt$accession == "MK153132"] <- "Arracacha virus B MK153132"
long.sim_nt$accession[long.sim_nt$accession == "MW582786"] <- "Cheravirus arracaciae MW582786"

long.sim_aa$accession[long.sim_aa$accession == "AB030941"] <- "Cheravirus mali AB030941"
long.sim_aa$accession[long.sim_aa$accession == "AJ621358"] <- "Cheravirus avii AJ621358"
long.sim_aa$accession[long.sim_aa$accession == "KT692953"] <- "Cheravirus ribis KT692953"
long.sim_aa$accession[long.sim_aa$accession == "DQ143875"] <- "Cheravirus pruni DQ143875"
long.sim_aa$accession[long.sim_aa$accession == "MK153132"] <- "Arracacha virus B MK153132"
long.sim_aa$accession[long.sim_aa$accession == "MW582786"] <- "Cheravirus arracaciae MW582786"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Cheravirus mali AB030941", "Cheravirus avii AJ621358",
                                                                  "Cheravirus ribis KT692953","Cheravirus pruni DQ143875",
                                                                  "Arracacha virus B MK153132","Cheravirus arracaciae MW582786"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Cheravirus mali AB030941", "Cheravirus avii AJ621358",
                                                                  "Cheravirus ribis KT692953","Cheravirus pruni DQ143875",
                                                                  "Arracacha virus B MK153132","Cheravirus arracaciae MW582786"))
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
cheravirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Pteropus rufus cheravirus RNA2 OQ818330")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

cheravirus_ictv_aa

## nucleotide
cheravirus_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

cheravirus_ictv_nt

library(cowplot)
cheravirus_ictv<-plot_grid(
  cheravirus_ictv_aa, cheravirus_ictv_nt,
  labels = NULL, ncol = 1
)
cheravirus_ictv



#Felisavirus
felisavirus_aa_ictv <- read.csv(file = "felisa_ictv_aa_partial_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
felisavirus_nt_ictv <- read.csv(file = "felisa_ictv_nt_partial_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(felisavirus_nt_ictv)
head(felisavirus_aa_ictv)

#move to long
long.sim_nt <- melt(felisavirus_nt_ictv, id.vars = c("pointer"), measure.vars = c("OQ818341","KX644943"))
long.sim_aa <- melt(felisavirus_aa_ictv, id.vars = c("pointer"), measure.vars = c("KX644943","OQ818341"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818341"] <- "Eidolon dupreanum felisavirus OQ818341"
long.sim_nt$accession[long.sim_nt$accession == "KX644943"] <- "Eidolon helvum felisavirus KX644943"
long.sim_aa$accession[long.sim_aa$accession == "OQ818341"] <- "Eidolon dupreanum felisavirus OQ818341"
long.sim_aa$accession[long.sim_aa$accession == "KX644943"] <- "Eidolon helvum felisavirus KX644943"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum felisavirus OQ818341", "Eidolon helvum felisavirus KX644943" 
))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon dupreanum felisavirus OQ818341", "Eidolon helvum felisavirus KX644943"
))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
felisavirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Pteropus rufus felisavirus OQ818335")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

felisavirus_ictv_aa

## nucleotide
felisavirus_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

felisavirus_ictv_nt

library(cowplot)
felisavirus_ictv<-plot_grid(
  felisavirus_ictv_aa, felisavirus_ictv_nt,
  labels = NULL, ncol = 1
)
felisavirus_ictv



#Hepatovirus
hepato_ictv_nt_partial <- read.csv(file = "hepato_ictv_nt_partial_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
hepato_ictv_aa_partial <- read.csv(file = "hepato_ictv_aa_partial_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(hepato_ictv_nt_partial)
head(hepato_ictv_aa_partial)

#move to long
long.sim_nt <- melt(hepato_ictv_nt_partial, id.vars = c("pointer"), measure.vars = c("HPA","KR703607","KT452637","KT452658",
                                                                                     "KT452685","KT452714","KT452730", "KT452735",
                                                                                     "KT452742"))
long.sim_aa <- melt(hepato_ictv_aa_partial, id.vars = c("pointer"), measure.vars = c("HPA","KR703607","KT452637","KT452658",
                                                                                     "KT452685","KT452714","KT452730", "KT452735",
                                                                                     "KT452742"))
unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "HPA"] <- "Hepatovirus A HPA"
long.sim_nt$accession[long.sim_nt$accession == "KR703607"] <- "Hepatovirus B KR703607"
long.sim_nt$accession[long.sim_nt$accession == "KT452637"] <- "Hepatovirus D KT452637"
long.sim_nt$accession[long.sim_nt$accession == "KT452658"] <- "Hepatovirus I KT452658"
long.sim_nt$accession[long.sim_nt$accession == "KT452685"] <- "Hepatovirus F KT452685"
long.sim_nt$accession[long.sim_nt$accession == "KT452714"] <- "Hepatovirus H KT452714"
long.sim_nt$accession[long.sim_nt$accession == "KT452730"] <- "Hepatovirus G KT452730"
long.sim_nt$accession[long.sim_nt$accession == "KT452735"] <- "Hepatovirus E KT452735"
long.sim_nt$accession[long.sim_nt$accession == "KT452742"] <- "Hepatovirus C KT452742"

long.sim_aa$accession[long.sim_aa$accession == "HPA"] <- "Hepatovirus A HPA"
long.sim_aa$accession[long.sim_aa$accession == "KR703607"] <- "Hepatovirus B KR703607"
long.sim_aa$accession[long.sim_aa$accession == "KT452637"] <- "Hepatovirus D KT452637"
long.sim_aa$accession[long.sim_aa$accession == "KT452658"] <- "Hepatovirus I KT452658"
long.sim_aa$accession[long.sim_aa$accession == "KT452685"] <- "Hepatovirus F KT452685"
long.sim_aa$accession[long.sim_aa$accession == "KT452714"] <- "Hepatovirus H KT452714"
long.sim_aa$accession[long.sim_aa$accession == "KT452730"] <- "Hepatovirus G KT452730"
long.sim_aa$accession[long.sim_aa$accession == "KT452735"] <- "Hepatovirus E KT452735"
long.sim_aa$accession[long.sim_aa$accession == "KT452742"] <- "Hepatovirus C KT452742"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Hepatovirus A HPA", "Hepatovirus B KR703607",
                                                                  "Hepatovirus C KT452742","Hepatovirus D KT452637",
                                                                  "Hepatovirus E KT452735","Hepatovirus F KT452685",
                                                                  "Hepatovirus G KT452730", "Hepatovirus H KT452714",
                                                                  "Hepatovirus I KT452658"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Hepatovirus A HPA", "Hepatovirus B KR703607",
                                                                  "Hepatovirus C KT452742","Hepatovirus D KT452637",
                                                                  "Hepatovirus E KT452735","Hepatovirus F KT452685",
                                                                  "Hepatovirus G KT452730", "Hepatovirus H KT452714",
                                                                  "Hepatovirus I KT452658"))
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
hepatovirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference:Eidolon dupreanum hepatovirus OQ818337")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_ictv_aa

## nucleotide
hepatovirus_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_ictv_nt

library(cowplot)
hepatovirus_ictv<-plot_grid(
  hepatovirus_ictv_aa, hepatovirus_ictv_nt,
  labels = NULL, ncol = 1
)
hepatovirus_ictv



#kobuvirus
kobuvirus_aa_ictv <- read.csv(file = "kobu_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
kobuvirus_nt_ictv <- read.csv(file = "kobu_ictv_nt_full_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(kobuvirus_nt_ictv)
head(kobuvirus_aa_ictv)

#move to long
long.sim_nt <- melt(kobuvirus_nt_ictv, id.vars = c("pointer"), measure.vars = c("AB040749","AB084788",
                                                                                "EU787450","KJ641686",
                                                                                "KT325853","LC055961",
                                                                                "OP287812"))
long.sim_aa <- melt(kobuvirus_aa_ictv, id.vars = c("pointer"), measure.vars = c("AB040749","AB084788",
                                                                                "EU787450","KJ641686",
                                                                                "KT325853","LC055961",
                                                                                "OP287812"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "AB040749"] <- "Aichivirus A AB040749"
long.sim_nt$accession[long.sim_nt$accession == "AB084788"] <- "Aichivirus B AB084788"
long.sim_nt$accession[long.sim_nt$accession == "EU787450"] <- "Aichivirus C EU787450"
long.sim_nt$accession[long.sim_nt$accession == "KJ641686"] <- "Aichivirus F KJ641686"
long.sim_nt$accession[long.sim_nt$accession == "KT325853"] <- "Aichivirus E KT325853"
long.sim_nt$accession[long.sim_nt$accession == "LC055961"] <- "Aichivirus D LC055961"
long.sim_nt$accession[long.sim_nt$accession == "OP287812"] <- "Eidolon dupreanum kobuvirus OP287812"

long.sim_aa$accession[long.sim_aa$accession == "AB040749"] <- "Aichivirus A AB040749"
long.sim_aa$accession[long.sim_aa$accession == "AB084788"] <- "Aichivirus B AB084788"
long.sim_aa$accession[long.sim_aa$accession == "EU787450"] <- "Aichivirus C EU787450"
long.sim_aa$accession[long.sim_aa$accession == "KJ641686"] <- "Aichivirus F KJ641686"
long.sim_aa$accession[long.sim_aa$accession == "KT325853"] <- "Aichivirus E KT325853"
long.sim_aa$accession[long.sim_aa$accession == "LC055961"] <- "Aichivirus D LC055961"
long.sim_aa$accession[long.sim_aa$accession == "OP287812"] <- "Eidolon dupreanum kobuvirus OP287812"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Aichivirus A AB040749", "Aichivirus B AB084788",
                                                                  "Aichivirus C EU787450","Aichivirus D LC055961",
                                                                  "Aichivirus E KT325853","Aichivirus F KJ641686",
                                                                  "Eidolon dupreanum kobuvirus OP287812"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Aichivirus A AB040749", "Aichivirus B AB084788",
                                                                  "Aichivirus C EU787450","Aichivirus D LC055961",
                                                                  "Aichivirus E KT325853","Aichivirus F KJ641686",
                                                                  "Eidolon dupreanum kobuvirus OP287812"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
kobuvirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Eidolon dupreanum kobuvirus OQ818322")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobuvirus_ictv_aa

## nucleotide
kobuvirus_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kobuvirus_ictv_nt

library(cowplot)
kobuvirus_ictv<-plot_grid(
  kobuvirus_ictv_aa, kobuvirus_ictv_nt,
  labels = NULL, ncol = 1
)
kobuvirus_ictv



#kunsagivirus
kunsagivirus_aa_ictv <- read.csv(file = "kun_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
kunsagivirus_nt_ictv <- read.csv(file = "kun_ictv_nt_full_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(kunsagivirus_nt_ictv)
head(kunsagivirus_aa_ictv)

#move to long
long.sim_nt <- melt(kunsagivirus_nt_ictv, id.vars = c("pointer"), measure.vars = c("KC935379","KX644936",
                                                                                "KY670597"))
long.sim_aa <- melt(kunsagivirus_aa_ictv, id.vars = c("pointer"), measure.vars = c("KC935379","KX644936",
                                                                                   "KY670597"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "KC935379"] <- "Kunsagivirus A KC935379"
long.sim_nt$accession[long.sim_nt$accession == "KX644936"] <- "Kunsagivirus B KX644936"
long.sim_nt$accession[long.sim_nt$accession == "KY670597"] <- "Kunsagivirus C KY670597"

long.sim_aa$accession[long.sim_aa$accession == "KC935379"] <- "Kunsagivirus A KC935379"
long.sim_aa$accession[long.sim_aa$accession == "KX644936"] <- "Kunsagivirus B KX644936"
long.sim_aa$accession[long.sim_aa$accession == "KY670597"] <- "Kunsagivirus C KY670597"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Kunsagivirus A KC935379", "Kunsagivirus B KX644936",
                                                                  "Kunsagivirus C KY670597"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Kunsagivirus A KC935379", "Kunsagivirus B KX644936",
                                                                  "Kunsagivirus C KY670597"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
kunsagivirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Eidolon dupreanum kunsagivirus OQ818317")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagivirus_ictv_aa

## nucleotide
kunsagivirus_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagivirus_ictv_nt

library(cowplot)
kunsagivirus_ictv<-plot_grid(
  kunsagivirus_ictv_aa, kunsagivirus_ictv_nt,
  labels = NULL, ncol = 1
)
kunsagivirus_ictv


#Mischivirus
mischivirus_aa_ictv <- read.csv(file = "mischi_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
mischivirus_nt_ictv <- read.csv(file = "mischi_ictv_nt_full_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(mischivirus_nt_ictv)
head(mischivirus_aa_ictv)

#move to long
long.sim_nt <- melt(mischivirus_nt_ictv, id.vars = c("pointer"), measure.vars = c("JQ814851","KP054273",
                                                                                  "KP100644","KY512802",
                                                                                  "MF352410"))
long.sim_aa <- melt(mischivirus_aa_ictv, id.vars = c("pointer"), measure.vars = c("JQ814851","KP054273",
                                                                                  "KP100644","KY512802",
                                                                                  "MF352410"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "JQ814851"] <- "Mischivirus A JQ814851"
long.sim_nt$accession[long.sim_nt$accession == "KP054273"] <- "Mischivirus B KP054273"
long.sim_nt$accession[long.sim_nt$accession == "KP100644"] <- "Mischivirus C KP100644"
long.sim_nt$accession[long.sim_nt$accession == "KY512802"] <- "Mischivirus D KY512802"
long.sim_nt$accession[long.sim_nt$accession == "MF352410"] <- "Mischivirus E MF352410"

long.sim_aa$accession[long.sim_aa$accession == "JQ814851"] <- "Mischivirus A JQ814851"
long.sim_aa$accession[long.sim_aa$accession == "KP054273"] <- "Mischivirus B KP054273"
long.sim_aa$accession[long.sim_aa$accession == "KP100644"] <- "Mischivirus C KP100644"
long.sim_aa$accession[long.sim_aa$accession == "KY512802"] <- "Mischivirus D KY512802"
long.sim_aa$accession[long.sim_aa$accession == "MF352410"] <- "Mischivirus E MF352410"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Mischivirus A JQ814851", "Mischivirus B KP054273",
                                                                  "Mischivirus C KP100644", "Mischivirus D KY512802",
                                                                  "Mischivirus E MF352410"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Mischivirus A JQ814851", "Mischivirus B KP054273",
                                                                  "Mischivirus C KP100644", "Mischivirus D KY512802",
                                                                  "Mischivirus E MF352410"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
mischivirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Pteropus rufus mischivirus OQ818316")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_ictv_aa

## nucleotide
mischivirus_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_ictv_nt

library(cowplot)
mischivirus_ictv<-plot_grid(
  mischivirus_ictv_aa, mischivirus_ictv_nt,
  labels = NULL, ncol = 1
)
mischivirus_ictv



#Sapelovirus full
sapelovirus_full_aa_ictv <- read.csv(file = "sapelo_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
sapelovirus_full_nt_ictv <- read.csv(file = "sapelo_ictv_nt_full_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(sapelovirus_full_nt_ictv)
head(sapelovirus_full_aa_ictv)

#move to long
long.sim_nt <- melt(sapelovirus_full_nt_ictv, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","AF406813","AY064708",
                                                                                  "NC_033820"))
long.sim_aa <- melt(sapelovirus_full_aa_ictv, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","AF406813","AY064708",
                                                                                       "NC_033820"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"


long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_nt$accession[long.sim_nt$accession == "AF406813"] <- "Sapelovirus A AF406813"
long.sim_nt$accession[long.sim_nt$accession == "AY064708"] <- "Sapelovirus B AY064708"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820"] <- "Eidolon helvum sapelovirus NC_033820"

long.sim_aa$accession[long.sim_aa$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_aa$accession[long.sim_aa$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_aa$accession[long.sim_aa$accession == "AF406813"] <- "Sapelovirus A AF406813"
long.sim_aa$accession[long.sim_aa$accession == "AY064708"] <- "Sapelovirus B AY064708"
long.sim_aa$accession[long.sim_aa$accession == "NC_033820"] <- "Eidolon helvum sapelovirus NC_033820"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321",
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Sapelovirus A AF406813", "Sapelovirus B AY064708",
                                                                  "Eidolon helvum sapelovirus NC_033820"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321",
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Sapelovirus A AF406813","Sapelovirus B AY064708",
                                                                  "Eidolon helvum sapelovirus NC_033820"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
sapelovirus_full_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Eidolon dupreanum sapelovirus OQ818320")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_full_ictv_aa

## nucleotide
sapelovirus_full_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_full_ictv_nt

library(cowplot)
sapelovirus_full_ictv<-plot_grid(
  sapelovirus_full_ictv_aa, sapelovirus_full_ictv_nt,
  labels = NULL, ncol = 1
)
sapelovirus_full_ictv


#Sapovirus partial
sapovirus_partial_aa_ictv <- read.csv(file = "sapo_ictv_aa_partial_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
sapovirus_partial_nt_ictv <- read.csv(file = "sapo_ictv_nt_partial_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(sapovirus_partial_nt_ictv)
head(sapovirus_partial_aa_ictv)

#move to long
long.sim_nt <- melt(sapovirus_partial_nt_ictv, id.vars = c("pointer"), measure.vars = c("OQ818340","OQ818347","OQ818348","HM002617","KX759623",
                                                                                       "OP963671"))
long.sim_aa <- melt(sapovirus_partial_aa_ictv, id.vars = c("pointer"), measure.vars = c("OQ818340","OQ818347","OQ818348","HM002617","KX759623",
                                                                                        "OP963671"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818340"] <- "Eidolon dupreanum sapovirus OQ818347"
long.sim_nt$accession[long.sim_nt$accession == "OQ818347"] <- "Rousettus madagascariensis sapovirus OQ818347"
long.sim_nt$accession[long.sim_nt$accession == "OQ818348"] <- "Rousettus madagascariensis sapovirus OQ818348"
long.sim_nt$accession[long.sim_nt$accession == "HM002617"] <- "Sapporo virus HM002617"
long.sim_nt$accession[long.sim_nt$accession == "KX759623"] <- "Bat sapovirus KX759623"
long.sim_nt$accession[long.sim_nt$accession == "OP963671"] <- "Rousettus leschenaultii sapelovirus OP963671"

long.sim_aa$accession[long.sim_aa$accession == "OQ818340"] <- "Eidolon dupreanum sapovirus OQ818347"
long.sim_aa$accession[long.sim_aa$accession == "OQ818347"] <- "Rousettus madagascariensis sapovirus OQ818347"
long.sim_aa$accession[long.sim_aa$accession == "OQ818348"] <- "Rousettus madagascariensis sapovirus OQ818348"
long.sim_aa$accession[long.sim_aa$accession == "HM002617"] <- "Sapporo virus HM002617"
long.sim_aa$accession[long.sim_aa$accession == "KX759623"] <- "Bat sapovirus KX759623"
long.sim_aa$accession[long.sim_aa$accession == "OP963671"] <- "Rousettus leschenaultii sapelovirus OP963671"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum sapovirus OQ818347",
                                                                  "Rousettus madagascariensis sapovirus OQ818347",
                                                                  "Rousettus madagascariensis sapovirus OQ818348",
                                                                  "Sapporo virus HM002617", "Bat sapovirus KX759623",
                                                                  "Rousettus leschenaultii sapelovirus OP963671"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Eidolon dupreanum sapovirus OQ818347",
                                                                  "Rousettus madagascariensis sapovirus OQ818347",
                                                                  "Rousettus madagascariensis sapovirus OQ818348",
                                                                  "Sapporo virus HM002617", "Bat sapovirus KX759623",
                                                                  "Rousettus leschenaultii sapelovirus OP963671"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
sapovirus_partial_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Eidolon dupreanum sapovirus OQ818319")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_partial_ictv_aa

## nucleotide
sapovirus_partial_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_partial_ictv_nt

library(cowplot)
sapovirus_partial_ictv<-plot_grid(
  sapovirus_partial_ictv_aa, sapovirus_partial_ictv_nt,
  labels = NULL, ncol = 1
)
sapovirus_partial_ictv



#Sapovirus full
sapovirus_full_aa_ictv <- read.csv(file = "sapo_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
sapovirus_full_nt_ictv <- read.csv(file = "sapo_ictv_nt_full_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(sapovirus_full_nt_ictv)
head(sapovirus_full_aa_ictv)

#move to long
long.sim_nt <- melt(sapovirus_full_nt_ictv, id.vars = c("pointer"), measure.vars = c("HM002617","KX759623"))
long.sim_aa <- melt(sapovirus_full_aa_ictv, id.vars = c("pointer"), measure.vars = c("HM002617","KX759623"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "HM002617"] <- "Sapporo virus HM002617"
long.sim_nt$accession[long.sim_nt$accession == "KX759623"] <- "Bat sapovirus KX759623"

long.sim_aa$accession[long.sim_aa$accession == "HM002617"] <- "Sapporo virus HM002617"
long.sim_aa$accession[long.sim_aa$accession == "KX759623"] <- "Bat sapovirus KX759623"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Sapporo virus HM002617", "Bat sapovirus KX759623"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Sapporo virus HM002617", "Bat sapovirus KX759623"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
sapovirus_full_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Eidolon dupreanum sapovirus OQ818319")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_full_ictv_aa

## nucleotide
sapovirus_full_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapovirus_full_ictv_nt

library(cowplot)
sapovirus_full_ictv<-plot_grid(
  sapovirus_full_ictv_aa, sapovirus_full_ictv_nt,
  labels = NULL, ncol = 1
)
sapovirus_full_ictv



#Teschovirus
teschovirus_aa_ictv <- read.csv(file = "tescho_ictv_aa_full_alignment.csv", header = T, stringsAsFactors = F) #nucleotide
teschovirus_nt_ictv <- read.csv(file = "tescho_ictv_nt_full_alignment.csv", header = T, stringsAsFactors = F) #animo acid
head(teschovirus_nt_ictv)
head(teschovirus_aa_ictv)

#move to long
long.sim_nt <- melt(teschovirus_nt_ictv, id.vars = c("pointer"), measure.vars = c("OQ818323","OQ818324","LC386158","MG875515","MT295502"))
long.sim_aa <- melt(teschovirus_aa_ictv, id.vars = c("pointer"), measure.vars = c("OQ818323","OQ818324","LC386158","MG875515","MT295502"))

unique(long.sim_nt$variable)
unique(long.sim_aa$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)
long.sim_aa$variable <- as.character(long.sim_aa$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"
names(long.sim_aa)[names(long.sim_aa)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818323"] <- "Rousettus madagascariensis teschovirus OQ818323"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324"] <- "Rousettus madagascariensis teschovirus OQ818324"
long.sim_nt$accession[long.sim_nt$accession == "LC386158"] <- "Teschovirus A LC386158"
long.sim_nt$accession[long.sim_nt$accession == "MG875515"] <- "Teschovirus B MG875515"
long.sim_nt$accession[long.sim_nt$accession == "MT295502"] <- "Teschovirus A6 MT295502"

long.sim_aa$accession[long.sim_aa$accession == "OQ818323"] <- "Rousettus madagascariensis teschovirus OQ818323"
long.sim_aa$accession[long.sim_aa$accession == "OQ818324"] <- "Rousettus madagascariensis teschovirus OQ818324"
long.sim_aa$accession[long.sim_aa$accession == "LC386158"] <- "Teschovirus A LC386158"
long.sim_aa$accession[long.sim_aa$accession == "MG875515"] <- "Teschovirus B MG875515"
long.sim_aa$accession[long.sim_aa$accession == "MT295502"] <- "Teschovirus A6 MT295502"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rousettus madagascariensis teschovirus OQ818323",
                                                                  "Rousettus madagascariensis teschovirus OQ818324",
                                                                  "Teschovirus A LC386158",
                                                                  "Teschovirus B MG875515",
                                                                  "Teschovirus A6 MT295502"))
long.sim_aa$accession <- factor(long.sim_aa$accession, levels = c("Rousettus madagascariensis teschovirus OQ818323",
                                                                  "Rousettus madagascariensis teschovirus OQ818324",
                                                                  "Teschovirus A LC386158",
                                                                  "Teschovirus B MG875515",
                                                                  "Teschovirus A6 MT295502"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <-0
long.sim_nt$value <- long.sim_nt$value/100

long.sim_aa$value[long.sim_aa$value<0] <- 0
long.sim_aa$value <- long.sim_aa$value/100

## amino acid
teschovirus_ictv_aa <- ggplot(long.sim_aa) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Amino acid similarity")+ xlab("") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = c(0.6,0.2), legend.direction = "vertical",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle("Reference: Eidolon dupreanum teschovirus OQ818318")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

teschovirus_ictv_aa

## nucleotide
teschovirus_ictv_nt <- ggplot(long.sim_nt) + geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("Nucleotide similarity")+xlab("Genome position") +
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position = "none", legend.direction = "horizontal",# legend.box = "vertical",
        legend.text = element_text(face="italic", size = 8),
        legend.title=element_text(face="italic", size = 8),
        axis.text = element_text(size=12), axis.title = element_text(size=12)) +
  #scale_color_manual(values=colz2)
  scale_fill_distiller()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

teschovirus_ictv_nt

library(cowplot)
teschovirus_ictv<-plot_grid(
  teschovirus_ictv_aa, teschovirus_ictv_nt,
  labels = NULL, ncol = 1
)
teschovirus_ictv
