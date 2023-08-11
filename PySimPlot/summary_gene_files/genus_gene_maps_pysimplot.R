library(ggplot2)
library(gggenes)
library(cowplot)
library(gridExtra)
library(grid)

homewd="/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/genome_annotation_and_characterization/"
setwd("~/Desktop/mada-bat-picornavirus/genome_annotation_and_characterization/genus_gene_maps")

#load files and make the genes into factor data
africa <- read.csv("africa_align_genes.csv", header = T, stringsAsFactors = F)
africa$gene<-factor(africa$gene, levels = c("5'UTR", "L peptide","VP4 peptide", "VP2 peptide", "VP0 peptide", "VP3 peptide",
                                               "VP1 peptide", "2A peptide", "2B peptide", "2C peptide", "3A peptide", "3B peptide",
                                               "3C peptide", "3D peptide", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein", 
                                               "Hypothetical protein", "Similar to structural polyprotein", 
                                              "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
blast <- read.csv("blast_only_genes.csv", header = T, stringsAsFactors = F)
blast$gene<-factor(blast$gene, levels = c("5'UTR", "L peptide","VP4 peptide", "VP2 peptide", "VP0 peptide", "VP3 peptide",
                                            "VP1 peptide", "2A peptide", "2B peptide", "2C peptide", "3A peptide", "3B peptide",
                                            "3C peptide", "3D peptide", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein", 
                                            "Hypothetical protein", "Similar to structural polyprotein", 
                                            "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
ictv <- read.csv("ictv_blast_genes.csv", header = T, stringsAsFactors = F)
ictv$gene<-factor(ictv$gene, levels = c("5'UTR", "L peptide","VP4 peptide", "VP2 peptide", "VP0 peptide", "VP3 peptide",
                                            "VP1 peptide", "2A peptide", "2B peptide", "2C peptide", "3A peptide", "3B peptide",
                                            "3C peptide", "3D peptide", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein", 
                                            "Hypothetical protein", "Similar to structural polyprotein", 
                                            "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))


##make the viruses into factor data
africa$molecule<-factor(africa$molecule, levels = c("Felisavirus", "Hepatovirus", "Kobuvirus", "Kunsagivirus",
                                                    "Mischivirus", "Sapovirus", "Sapelovirus"))
blast$molecule<-factor(blast$molecule, levels = c("Bat picornavirus","Cheravirus", "Felisavirus", "Hepatovirus", "Kobuvirus", "Kunsagivirus",
                                                    "Mischivirus", "Sapovirus", "Sapelovirus", "Teschovirus"))
ictv$molecule<-factor(ictv$molecule, levels = c("Bat picornavirus","Cheravirus", "Felisavirus", "Hepatovirus", "Kobuvirus", "Kunsagivirus",
                                                  "Mischivirus", "Sapovirus", "Sapelovirus", "Teschovirus"))



#Pick colors for genes
colz=c("5'UTR"="gold", "L peptide"="royalblue","VP4 peptide"="paleturquoise3", "VP2 peptide"="skyblue1", "VP0 peptide"="royalblue4", "VP3 peptide"="steelblue1",
        "VP1 peptide"="cadetblue1", "2A peptide"="palevioletred1", "2B peptide"="red4", "2C peptide"="palevioletred3", "3A peptide"="tomato2", "3B peptide"="plum",
        "3C peptide"="rosybrown1", "3D peptide"="pink2", 
        "Polyprotein"="azure3","Putative polyprotein"="mediumorchid1", "Non-structural polyprotein"="mediumorchid4", 
        "Putative minor structural protein" ="slateblue3",
        "Hypothetical protein"="darkslategrey", "Similar to structural polyprotein"="mediumpurple1", "Similar to putative polyprotein"="mediumpurple3", 
        "Similar to polyprotein"="mediumpurple4","3'UTR"="yellow")














