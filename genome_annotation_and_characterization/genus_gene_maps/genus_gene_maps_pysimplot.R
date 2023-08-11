rm(list=ls())

library(ggplot2)
library(gggenes)
library(cowplot)
library(gridExtra)
library(grid)

homewd="/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/genome_annotation_and_characterization/"
setwd("~/Desktop/mada-bat-picornavirus/genome_annotation_and_characterization/genus_gene_maps")

#load files and make the genes into factor data
africa <- read.csv("africa_align_genes.csv", header = T, stringsAsFactors = F)
africa$gene<-factor(africa$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                               "VP1", "2A", "2B", "2C", "3A", "3B",
                                               "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                               "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                              "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
blast <- read.csv("blast_only_genes.csv", header = T, stringsAsFactors = F)
blast$gene<-factor(blast$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                          "VP1", "2A", "2B", "2C", "3A", "3B",
                                          "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                          "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                          "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
ictv <- read.csv("ictv_blast_genes.csv", header = T, stringsAsFactors = F)
ictv$gene<-factor(ictv$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
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
blast_pep <- read.csv("blast_only_peptides.csv", header = T, stringsAsFactors = F)
blast_pep$gene<-factor(blast_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                  "VP1", "2A", "2B", "2C", "3A", "3B",
                                                  "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                                  "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                  "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
ictv_pep <- read.csv("ictv_blast_peptides.csv", header = T, stringsAsFactors = F)
ictv_pep$gene<-factor(ictv_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
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
blast_feat <- read.csv("blast_only_features.csv", header = T, stringsAsFactors = F)
blast_feat$gene<-factor(blast_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                  "VP1", "2A", "2B", "2C", "3A", "3B",
                                                  "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                                  "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                  "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
ictv_feat <- read.csv("ictv_blast_features.csv", header = T, stringsAsFactors = F)
ictv_feat$gene<-factor(ictv_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
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
  # geom_feature(data=africa_felisavirus_feat, aes(x=from, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=africa_felisavirus_feat, aes(x=from, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,9188)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_felis
  

africa_hep<-ggplot(africa_hepatovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=africa_hepatovirus_feat, aes(x=from, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=africa_hepatovirus_feat, aes(x=from, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7986)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_hep


africa_kobu<-ggplot(africa_kobuvirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=africa_kobuvirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=africa_kobuvirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = africa_kobuvirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = africa_kobuvirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8379)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_kobu


africa_kun<-ggplot(africa_kunsagivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=africa_kunsagivirus_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=africa_kunsagivirus_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7333)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_kun


africa_mischi<-ggplot(africa_mischivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=africa_mischivirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=africa_mischivirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = africa_mischivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = africa_mischivirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8886)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_mischi


africa_sapo<-ggplot(africa_sapovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=africa_sapovirus_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=africa_sapovirus_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7759)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_sapo


africa_sapelo<-ggplot(africa_sapelovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=africa_sapelovirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=africa_sapelovirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = africa_sapelovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = africa_sapelovirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7940)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
africa_sapelo








#plot BLAST only gene maps

#Subset BLAST data by virus
blast_batpicornavirus<-subset(blast,molecule=="Bat picornavirus")
blast_cheravirus<-subset(blast,molecule=="Cheravirus")
blast_felisavirus<-subset(blast,molecule=="Felisavirus")
blast_hepatovirus<-subset(blast,molecule=="Hepatovirus")
blast_kobuvirus<-subset(blast,molecule=="Kobuvirus")
blast_kunsagivirus<-subset(blast,molecule=="Kunsagivirus")
blast_mischivirus<-subset(blast,molecule=="Mischivirus")
blast_sapovirus<-subset(blast,molecule=="Sapovirus")
blast_sapovirus_full<-subset(blast_sapovirus,type=="full")
blast_sapovirus_partial<-subset(blast_sapovirus,type=="partial")
blast_sapelovirus<-subset(blast,molecule=="Sapelovirus")
blast_sapelovirus_full<-subset(blast_sapelovirus,type=="full")
blast_sapelovirus_partial<-subset(blast_sapelovirus,type=="partial")
blast_teschovirus<-subset(blast,molecule=="Teschovirus")

blast_batpicornavirus_pep<-subset(blast_pep,molecule=="Bat picornavirus")
blast_kobuvirus_pep<-subset(blast_pep,molecule=="Kobuvirus")
blast_mischivirus_pep<-subset(blast_pep,molecule=="Mischivirus")
blast_sapelovirus_pep<-subset(blast_pep,molecule=="Sapelovirus")
blast_sapelovirus_pep_full<-subset(blast_sapelovirus_pep,type=="full")
blast_sapelovirus_pep_partial<-subset(blast_sapelovirus_pep,type=="partial")
blast_teschovirus_pep<-subset(blast_pep,molecule=="Teschovirus")

blast_batpicornavirus_feat<-subset(blast_feat,molecule=="Bat picornavirus")
blast_cheravirus_feat<-subset(blast_feat,molecule=="Cheravirus")
blast_felisavirus_feat<-subset(blast_feat,molecule=="Felisavirus")
blast_hepatovirus_feat<-subset(blast_feat,molecule=="Hepatovirus")
blast_kobuvirus_feat<-subset(blast_feat,molecule=="Kobuvirus")
blast_kunsagivirus_feat<-subset(blast_feat,molecule=="Kunsagivirus")
blast_mischivirus_feat<-subset(blast_feat,molecule=="Mischivirus")
blast_sapovirus_feat<-subset(blast_feat,molecule=="Sapovirus")
blast_sapovirus_full_feat<-subset(blast_sapovirus_feat,type=="full")
blast_sapovirus_partial_feat<-subset(blast_sapovirus_feat,type=="partial")
blast_sapelovirus_feat<-subset(blast_feat,molecule=="Sapelovirus")
blast_sapelovirus_full_feat<-subset(blast_sapelovirus_feat,type=="full")
blast_sapelovirus_partial_feat<-subset(blast_sapelovirus_feat,type=="partial")
blast_teschovirus_feat<-subset(blast_feat,molecule=="Teschovirus")



#plot blast gene maps
blast_batpicorna<-ggplot(blast_batpicornavirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=blast_batpicornavirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=blast_batpicornavirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = blast_batpicornavirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = blast_batpicornavirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8079)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_batpicorna


blast_chera<-ggplot(blast_cheravirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=blast_cheravirus_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=blast_cheravirus_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,3510)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_chera


blast_felisa<-ggplot(blast_felisavirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=blast_felisavirus_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=blast_felisavirus_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,9188)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_felisa


blast_hepato<-ggplot(blast_hepatovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=blast_hepatovirus_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=blast_hepatovirus_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7485)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_hepato


blast_kobu<-ggplot(blast_kobuvirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=blast_kobuvirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=blast_kobuvirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = blast_kobuvirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                   xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = blast_kobuvirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8379)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_kobu


blast_kun<-ggplot(blast_kunsagivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=blast_kunsagivirus_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=blast_kunsagivirus_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7333)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_kun


blast_mischi<-ggplot(blast_mischivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=blast_mischivirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=blast_mischivirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = blast_mischivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                             xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = blast_mischivirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8588)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_mischi


blast_sapelo_full<-ggplot(blast_sapelovirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=blast_sapelovirus_full, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=blast_sapelovirus_full, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = blast_sapelovirus_pep_full, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                               xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = blast_sapelovirus_pep_full, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7936)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_sapelo_full


blast_sapelo_partial<-ggplot(blast_sapelovirus_partial, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=blast_sapelovirus_partial_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=blast_sapelovirus_partial_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7936)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_sapelo_partial


blast_sapo_full<-ggplot(blast_sapovirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=blast_sapovirus_full_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=blast_sapovirus_full_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7573)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_sapo_full


blast_sapo_partial<-ggplot(blast_sapovirus_partial, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=blast_sapovirus_full_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=blast_sapovirus_full_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7573)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_sapo_partial


blast_tescho<-ggplot(blast_teschovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=blast_teschovirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=blast_teschovirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = blast_teschovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                    xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = blast_teschovirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7695)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
blast_tescho







#Plot ICTV and BLAST together plots

#Subset ICTV data by virus
ictv_batpicornavirus<-subset(ictv,molecule=="Bat picornavirus")
ictv_cheravirus<-subset(ictv,molecule=="Cheravirus")
ictv_felisavirus<-subset(ictv,molecule=="Felisavirus")
ictv_hepatovirus<-subset(ictv,molecule=="Hepatovirus")
ictv_kobuvirus<-subset(ictv,molecule=="Kobuvirus")
ictv_kunsagivirus<-subset(ictv,molecule=="Kunsagivirus")
ictv_mischivirus<-subset(ictv,molecule=="Mischivirus")
ictv_sapovirus<-subset(ictv,molecule=="Sapovirus")
ictv_sapovirus_full<-subset(ictv_sapovirus,type=="full")
ictv_sapovirus_partial<-subset(ictv_sapovirus,type=="partial")
ictv_sapelovirus<-subset(ictv,molecule=="Sapelovirus")
ictv_sapelovirus_full<-subset(ictv_sapelovirus,type=="full")
ictv_sapelovirus_partial<-subset(ictv_sapelovirus,type=="partial")
ictv_teschovirus<-subset(ictv,molecule=="Teschovirus")

ictv_batpicornavirus_pep<-subset(ictv_pep,molecule=="Bat picornavirus")
ictv_kobuvirus_pep<-subset(ictv_pep,molecule=="Kobuvirus")
ictv_mischivirus_pep<-subset(ictv_pep,molecule=="Mischivirus")
ictv_sapelovirus_pep<-subset(ictv_pep,molecule=="Sapelovirus")
ictv_sapelovirus_pep_full<-subset(ictv_sapelovirus_pep,type=="full")
ictv_sapelovirus_pep_partial<-subset(ictv_sapelovirus_pep,type=="partial")
ictv_teschovirus_pep<-subset(ictv_pep,molecule=="Teschovirus")

ictv_batpicornavirus_feat<-subset(ictv_feat,molecule=="Bat picornavirus")
ictv_cheravirus_feat<-subset(ictv_feat,molecule=="Cheravirus")
ictv_felisavirus_feat<-subset(ictv_feat,molecule=="Felisavirus")
ictv_hepatovirus_feat<-subset(ictv_feat,molecule=="Hepatovirus")
ictv_kobuvirus_feat<-subset(ictv_feat,molecule=="Kobuvirus")
ictv_kunsagivirus_feat<-subset(ictv_feat,molecule=="Kunsagivirus")
ictv_mischivirus_feat<-subset(ictv_feat,molecule=="Mischivirus")
ictv_sapovirus_feat<-subset(ictv_feat,molecule=="Sapovirus")
ictv_sapovirus_full_feat<-subset(ictv_sapovirus_feat,type=="full")
ictv_sapovirus_partial_feat<-subset(ictv_sapovirus_feat,type=="partial")
ictv_sapelovirus_feat<-subset(ictv_feat,molecule=="Sapelovirus")
ictv_sapelovirus_full_feat<-subset(ictv_sapelovirus_feat,type=="full")
ictv_sapelovirus_partial_feat<-subset(ictv_sapelovirus_feat,type=="partial")
ictv_teschovirus_feat<-subset(ictv_feat,molecule=="Teschovirus")


#plot ictv and blast plots
ictv_batpicorna<-ggplot(ictv_batpicornavirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=ictv_batpicornavirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=ictv_batpicornavirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_batpicornavirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                   xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = ictv_batpicornavirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8331)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_batpicorna


ictv_chera<-ggplot(ictv_cheravirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_cheravirus_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_cheravirus_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,3843)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_chera


ictv_felisa<-ggplot(ictv_felisavirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_felisavirus_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_felisavirus_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,9188)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_felisa


ictv_hepato<-ggplot(ictv_hepatovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8441)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_hepato


ictv_kobu<-ggplot(ictv_kobuvirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=ictv_kobuvirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=ictv_kobuvirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_kobuvirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                             xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = ictv_kobuvirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,9377)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_kobu


ictv_kun<-ggplot(ictv_kunsagivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_kunsagivirus_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_kunsagivirus_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8098)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_kun


ictv_mischi<-ggplot(ictv_mischivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=ictv_mischivirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=ictv_mischivirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_mischivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                               xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = ictv_mischivirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,9572)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_mischi


ictv_sapelo_full<-ggplot(ictv_sapelovirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=ictv_sapelovirus_full, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=ictv_sapelovirus_full, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_sapelovirus_pep_full, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                    xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = ictv_sapelovirus_pep_full, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8650)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapelo_full


ictv_sapelo_partial<-ggplot(ictv_sapelovirus_partial, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_sapelovirus_partial_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_sapelovirus_partial_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(8650)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapelo_partial


ictv_sapo_full<-ggplot(ictv_sapovirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_sapovirus_full_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_sapovirus_full_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8017)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapo_full


ictv_sapo_partial<-ggplot(ictv_sapovirus_partial, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_sapovirus_full_feat, aes(x=start, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_sapovirus_full_feat, aes(x=start, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_gene_label(aes(label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,8037)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapo_partial


ictv_tescho<-ggplot(ictv_teschovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  geom_feature(data=ictv_teschovirus, aes(x=start, y=molecule),
               feature_height = grid::unit(6,"mm"))+
  geom_feature_label(data=ictv_teschovirus, aes(x=start, y=molecule, label=gene),
                     feature_height = grid::unit(6,"mm"),
                     label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_teschovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                               xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_subgene_label(data = ictv_teschovirus_pep, aes(xmin = from, xmax = to,xsubmin=from, xsubmax=to, label=gene), align = "left")+
  scale_fill_manual(values=colz)+
  theme_genes()+
  xlim(0,7725)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_tescho

