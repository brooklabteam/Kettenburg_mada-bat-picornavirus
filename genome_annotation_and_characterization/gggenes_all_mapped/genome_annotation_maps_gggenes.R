library(ggplot2)
library(gggenes)
library(cowplot)
library(gridExtra)
library(grid)

homewd="/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/genome_annotation_and_characterization/"


setwd("~/Desktop/mada-bat-picornavirus/genome_annotation_and_characterization")
##Without polyprotein, get the base information loaded
full_genome_gggenes_nopoly <- read.csv("full_genome_nopoly.csv", header = T, stringsAsFactors = F)
full_genome_gggenes_nopoly$gene<-factor(full_genome_gggenes_nopoly$gene, levels = c("5'UTR", "L peptide","VP4 peptide", "VP2 peptide", "VP0 peptide", "VP3 peptide",
                                                                                        "VP1 peptide", "2A peptide", "2B peptide", "2C peptide", "3A peptide", "3B peptide",
                                                                                        "3C peptide", "3D peptide", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein", 
                                                                                        "Hypothetical protein", "Similar to structural polyprotein", 
                                                                                        "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))

##With polyprotein, add on the peptides on top of the polyprotein figures
full_genome_gggenes_withpoly <- read.csv("full_genome_poly.csv", header = T, stringsAsFactors = F)

# colz=c("5'UTR"="gold", "L peptide"="steelblue1","VP4 peptide"="thistle1", "VP2 peptide"="firebrick1", "VP0 peptide"="sienna3", "VP3 peptide"="tomato1",
#        "VP1 peptide"="salmon", "2A peptide"="palevioletred1", "2B peptide"="red4", "2C peptide"="indianred1", "3A peptide"="plum1", "3B peptide"="lightskyblue1",
#        "3C peptide"="rosybrown1", "3D peptide"="skyblue3", 
#        "Polyprotein"="white","Putative polyprotein"="cadetblue2", "Non-structural polyprotein"="cadetblue4", 
#        "Putative minor structural protein" ="paleturquoise1",
#        "Hypothetical protein"="darkslategrey", "Similar to structural polyprotein"="slateblue1", "Similar to putative polyprotein"="slateblue3", 
#        "Similar to polyprotein"="slateblue4","3'UTR"="yellow")

colz2=c("5'UTR"="gold", "L peptide"="royalblue","VP4 peptide"="paleturquoise3", "VP2 peptide"="skyblue1", "VP0 peptide"="royalblue4", "VP3 peptide"="steelblue1",
       "VP1 peptide"="cadetblue1", "2A peptide"="palevioletred1", "2B peptide"="red4", "2C peptide"="palevioletred3", "3A peptide"="tomato2", "3B peptide"="plum",
       "3C peptide"="rosybrown1", "3D peptide"="pink2", 
       "Polyprotein"="azure3","Putative polyprotein"="mediumorchid1", "Non-structural polyprotein"="mediumorchid4", 
       "Putative minor structural protein" ="slateblue3",
       "Hypothetical protein"="darkslategrey", "Similar to structural polyprotein"="mediumpurple1", "Similar to putative polyprotein"="mediumpurple3", 
       "Similar to polyprotein"="mediumpurple4","3'UTR"="yellow")

full_genome_gggenes_withpoly$gene<-factor(full_genome_gggenes_withpoly$gene, levels = c("5'UTR", "L peptide","VP4 peptide", "VP2 peptide", "VP0 peptide", "VP3 peptide",
                                                                                        "VP1 peptide", "2A peptide", "2B peptide", "2C peptide", "3A peptide", "3B peptide",
                                                                                        "3C peptide", "3D peptide", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein", 
                                                                                        "Hypothetical protein", "Similar to structural polyprotein", 
                                                                                        "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))

full_genome_gggenes_withpoly$virus<-factor(full_genome_gggenes_withpoly$virus, levels=c("Mischivirus", "Kobuvirus", "Kunsagivirus", "Nepovirus", "Picornavirus", "Sapelovirus",
                                                                                        "Sapovirus", "Teschovirus"))


full1<-ggplot(full_genome_gggenes_withpoly, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(2, "mm"),
                  arrowhead_height = grid::unit(2, "mm"),
                  arrow_body_height = grid::unit(2, "mm")) +
  facet_wrap(~ virus, scales = "free", ncol = 1) +
  geom_subgene_arrow(full_genome_gggenes_nopoly, mapping=aes(xmin = start, xmax = end, y = molecule, fill=gene,
                                                     xsubmin=start, xsubmax=end), color="black", alpha=.7,
                     arrowhead_width = grid::unit(2, "mm"),
                     arrowhead_height = grid::unit(2, "mm"),
                     arrow_body_height = grid::unit(2, "mm"))+
  scale_fill_manual(values=colz2)+
  theme_linedraw()+
  theme(legend.position = "none")+
  xlim(0,8500)+
  xlab("Genome position") + ylab("Accession")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

full1
#7x5 inch PDF portrait 1 column
#7x4 inch PDF landscape 2 column








##Partial genomes
partial_genome_gggenes <- read.csv("partial_genome.csv", header = T, stringsAsFactors = F)
partial_genome_gggenes$gene<-factor(partial_genome_gggenes$gene, levels = c("5'UTR", "L peptide","VP4 peptide", "VP2 peptide", "VP0 peptide", "VP3 peptide",
                                                                            "VP1 peptide", "2A peptide", "2B peptide", "2C peptide", "3A peptide", "3B peptide",
                                                                            "3C peptide", "3D peptide", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein", 
                                                                            "Hypothetical protein", "Similar to structural polyprotein", 
                                                                            "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))

partial_genome_gggenes$virus<-factor(partial_genome_gggenes$virus, levels=c("Bee virus", "Cheravirus", "Cripavirus", "Felisavirus", "Hepatovirus", "Picorna-like virus", "Picornavirus", "Sapelovirus",
                                                                            "Sapovirus", "Tetnovirus"))


partial1<-ggplot(partial_genome_gggenes, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(2, "mm"),
                  arrowhead_height = grid::unit(2, "mm"),
                  arrow_body_height = grid::unit(2, "mm")) +
  facet_wrap(~ virus, scales = "free", ncol = 1) +
  scale_fill_manual(values=colz2)+
  theme_linedraw()+
  theme(legend.position = "none")+
  xlim(0,8500)+
  xlab("Genome position") + ylab("Accession")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
partial1
#7x4 landscape 2 column
#10x5 inch PDF portrait 1 column






##all genes on one legend
setwd("~/Desktop/mada-bat-picornavirus/genome_annotation_and_characterization")
legend<-read.csv("all_for_legend.csv", header = T, stringsAsFactors = F)
legend$gene<-factor(legend$gene, levels = c("5'UTR", "L peptide","VP4 peptide", "VP2 peptide", "VP0 peptide", "VP3 peptide",
                                            "VP1 peptide", "2A peptide", "2B peptide", "2C peptide", "3A peptide", "3B peptide",
                                            "3C peptide", "3D peptide", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein", 
                                            "Hypothetical protein", "Similar to structural polyprotein", 
                                            "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))

legend$virus<-factor(legend$virus, levels=c("Mischivirus", 
                                            "Kobuvirus", 
                                            "Kunsagivirus", 
                                            "Nepovirus", 
                                            "Picornavirus", 
                                            "Sapelovirus",
                                            "Sapovirus", 
                                            "Teschovirus",
                                            "Bee virus", 
                                            "Cheravirus", 
                                            "Cripavirus", 
                                            "Felisavirus", 
                                            "Hepatovirus", 
                                            "Picorna-like virus",
                                            "Tetnovirus"))


legend_fig<-ggplot(legend, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(3, "mm"),
                  arrow_body_height = grid::unit(3, "mm")) +
  facet_wrap(~ virus, scales = "free", ncol = 3) +
  scale_fill_manual(values=colz2)+
  theme_linedraw()+
  theme(legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.position = "right")+
  xlim(0,8500)+
  xlab("Genome position") + ylab("Accession")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
legend_fig

#get the legend separately
#install.packages('ggpubr')
library(ggpubr)
leg_bottom<-get_legend(legend_fig)
leg_bottom<-as_ggplot(leg_bottom)
leg_bottom <- leg_bottom + theme(
  legend.margin=margin(c(0,0,0,0)))

leg_right<-get_legend(legend_fig)
leg_right<-as_ggplot(leg_right)
leg_right <- leg_right + theme(
  legend.margin=margin(c(0,0,0,0)))
#4x7 inch PDF landscape bottom legend
#3x3.1 inch PDF portrait right legend

##Can't get the legend to go nicely into plot_grid so I'll add the legend later separately

#Put it all together
seq1<-cowplot::plot_grid(full1,partial1,ncol=2, labels = "AUTO")
seq1

seq2<-cowplot::plot_grid(full2,partial2,ncol=2, labels = "AUTO")
seq2
