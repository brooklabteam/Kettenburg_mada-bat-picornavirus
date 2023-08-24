#Part A phylogenetic tree

rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)
library(phylotools)
library(phylobase)
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

homewd= "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/"

##Picornavirales all full and partial sequences >3kb
setwd(paste0(homewd,"/raxml_trees/bat_picornavirales/all"))

#load the tree
tree <-  read.tree("T10.raxml.supportFBP") 

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("bat_picornavirales_all_over_3kb.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #256
length(tree$tip.label) #256

#check subgroup names
unique(dat$Family)

colz = c("Caliciviridae" = "royalblue3",    "Picornaviridae"  = "turquoise1",   "Secoviridae"  = "goldenrod1",
         "Coronaviridae"  = "black", "Unclassified"  = "darkorange1")

#pick order for the labels
dat$Family <- factor(dat$Family, levels = c("Secoviridae","Picornaviridae",
                                            "Caliciviridae", "Unclassified","Coronaviridae")) 

dat$novel <- as.factor(dat$novel)


#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Family)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now check that you don't have NAs and blanks throughout the component parts of the name
#you also can edit the original csv file to replace these blanks if you can find the correct info

dat$Isolate  #some are blank, so convert those to NA
dat$Isolate[dat$Isolate==""] <- NA
dat$Accession #all good

dat$source
dat$source[dat$source==""] <- NA
dat$Country #some are blank, so convert those to NA
dat$Country[dat$Country==""] <- NA
dat$Collection_Date #these are messy, some are years and some are full dates. I just want years, will manually fix


#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Genus[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and source only:
dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Genus[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Genus[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                                             #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])
#look at dat$new_label
dat$new_label #they all look great

#make sure to sort in order
colnames(dat)[colnames(dat)=="Accession"]="old_tip_label"
tree.dat <- data.frame(old_tip_label=rooted.tree$tip.label, num =1:length(rooted.tree$tip.label))
head(tree.dat)
head(dat)
tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = F, sort = F)

names(tree.dat)

tree.dat$tip_label <- tree.dat$new_label
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Family, Genus, novel, old_tip_label)

rooted.tree$tip.label <- tree.dat$tip_label

head(tree.dat)
#verify that the old labels match the new

cbind(tree.dat$old_tip_label, rooted.tree$tip.label)

#check out the labels
tree.dat$tip_label#all look good

tree.dat$Host[tree.dat$Host==0] <- "Non-bat host"
tree.dat$Host[tree.dat$Host==1] <- "Bat host"
tree.dat$Host <- as.factor(tree.dat$Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Family, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none", shape="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1.5, linesize = .5) +
  theme(legend.margin = margin(),
        legend.position = "left", 
        legend.direction = "vertical",
        legend.text = element_text(size=8), 
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  xlim(c(0,20))

p1

##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree
#collapse the labeled clades
p3<-collapse(p1, 353)+geom_point2(aes(subset=(node==353)), size=3, shape=22, fill="turquoise1")
p4<-collapse(p3, 384)+geom_point2(aes(subset=(node==384)), size=3, shape=22, fill="turquoise1")
p5<-collapse(p4, 407)+geom_point2(aes(subset=(node==407)), size=3, shape=22, fill="turquoise1")
p6<-collapse(p5, 422)+geom_point2(aes(subset=(node==422)), size=3, shape=22, fill="turquoise1")
p7<-collapse(p6, 315)+geom_point2(aes(subset=(node==315)), size=3, shape=22, fill="turquoise1")
p8<-collapse(p7, 461)+geom_point2(aes(subset=(node==461)), size=3, shape=22, fill="turquoise1")
p9<-collapse(p8, 281)+geom_point2(aes(subset=(node==281)), size=3, shape=22, fill="turquoise1")
p9

p9<-p9+geom_cladelabel(node = 353, label = "Parechovirus/Scotophilus_kuhlii/Vietnam clade",offset=0.1, fontsize = 3, color="black") +
  geom_cladelabel(node = 384, label = "Shanbavirus/Rhinolophus_sp./China clade", offset=0.1,fontsize = 3, color="black") +
  geom_cladelabel(node = 407, label = "Shanbavirus/Miniopterus_shreibersii/Hungary clade", offset=0.1, fontsize = 3, color="black") +
  geom_cladelabel(node = 422, label = "Kobuvirus/Scotophilus_kuhlii/Vietnam clade", offset=0.1,fontsize = 3, color="black") +
  geom_cladelabel(node = 315, label = "Unclassified/Miniopterus_pusillus/Chian clade", offset=0.1,fontsize = 3, color="black") +
  geom_cladelabel(node = 461, label = "Shanbavirus/Rhinolophus_sp./China clade", offset=0.1, fontsize = 3, color="black") +
  geom_cladelabel(node = 281, label = "Unclassified/Scotophilus_kuhlii/Vietnam clade", offset=0.1, fontsize = 3, color="black")
p9


##add bootstrap values to this tree
p9.dat <- p9$data
p9.dat$Bootstrap <- NA
Bootstrap<-p9.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p9.dat$label)] <- as.numeric(p9.dat$label[(length(tree.dat$tip_label)+1):length(p9.dat$label)])#fill with label

batpicornavirales_tree <- p9  %<+% p9.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
batpicornavirales_tree





 #Now plot the bootscan and gene maps for part B
##First make the genome plots
homewd="/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/genome_annotation_and_characterization/"
setwd("~/Desktop/developer/mada-bat-picornavirus/genome_annotation_and_characterization/genus_gene_maps")

#Load the gene data
bat <- read.csv("bat_genes.csv", header = T, stringsAsFactors = F)
bat$gene<-factor(bat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                      "VP1", "2A", "2B", "2C", "3A", "3B",
                                      "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                      "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                      "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the peptide files
bat_pep <- read.csv("bat_peptides.csv", header = T, stringsAsFactors = F)
bat_pep$gene<-factor(bat_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                              "VP1", "2A", "2B", "2C", "3A", "3B",
                                              "3C", "3D", "Polyprotein", "Putative polyprotein", "Putative minor structural protein", "Non-structural polyprotein",
                                              "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                              "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the feature file in case its needed
bat_feat <- read.csv("bat_features.csv", header = T, stringsAsFactors = F)
bat_feat$gene<-factor(bat_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
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


#Plot gene maps first

#Subset bat data by virus and dataset (all or full bat picornavirales)
bat_hepatovirus<-subset(bat,molecule=="Hepatovirus")
bat_hepatovirus_all<-subset(bat_hepatovirus,type=="all")
bat_mischivirus<-subset(bat,molecule=="Mischivirus")
bat_mischivirus_all<-subset(bat_mischivirus,type=="all")
bat_sapelovirus<-subset(bat,molecule=="Sapelovirus")
bat_sapelovirus_all<-subset(bat_sapelovirus,type=="all")
bat_felispicornalike<-subset(bat,molecule=="Felisavirus and Picorna-like virus")

bat_hepatovirus_pep<-subset(bat_pep,molecule=="Hepatovirus")
bat_hepatovirus_all_pep<-subset(bat_hepatovirus_pep,type=="all")
bat_mischivirus_pep<-subset(bat_pep,molecule=="Mischivirus")
bat_mischivirus_all_pep<-subset(bat_mischivirus_pep,type=="all")
bat_sapelovirus_pep<-subset(bat_pep,molecule=="Sapelovirus")
bat_sapelovirus_all_pep<-subset(bat_sapelovirus_pep,type=="all")

bat_hepatovirus_feat<-subset(bat_feat,molecule=="Hepatovirus")
bat_hepatovirus_all_feat<-subset(bat_hepatovirus_feat,type=="all")
bat_mischivirus_feat<-subset(bat_feat,molecule=="Mischivirus")
bat_mischivirus_all_feat<-subset(bat_mischivirus_feat,type=="all")
bat_sapelovirus_feat<-subset(bat_feat,molecule=="Sapelovirus")
bat_sapelovirus_all_feat<-subset(bat_sapelovirus_feat,type=="all")
bat_felispicornalike_feat<-subset(bat_feat,molecule=="Felisavirus and Picorna-like virus")



##Now the bat gene maps



#all seq over 3kb

colzpalette<-c("#8ECAE6","#219EBC","#023047","#FFB703","#FB8500","#E48B97","#B52B09","#A60067","#987B6F","#8FD694")

felispicorna_all<-ggplot(bat_felispicornalike, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_text(data=bat_felispicornalike_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(3600,11970),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
felispicorna_all


bat_hepato_all<-ggplot(bat_hepatovirus_all, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_hepatovirus_all_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_hepatovirus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(800,6362),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_hepato_all


bat_mischi_all<-ggplot(bat_mischivirus_all, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_mischivirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_mischivirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_mischivirus_all_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_mischivirus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1590,8850),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_mischi_all


bat_sapelo_all<-ggplot(bat_sapelovirus_all, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=bat_sapelovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=bat_sapelovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = bat_sapelovirus_all_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=bat_sapelovirus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1290,8200),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
bat_sapelo_all

##Now get the plots from bootscan

setwd("~/Desktop/developer/mada-bat-picornavirus/bootscan/output_bat/all")

#Felisavirus
felisa_picornalike_bat_boot <- read.csv(file = "felisa_picorna_like_boot.csv", header = T, stringsAsFactors = F) 
head(felisa_picornalike_bat_boot)

#move to long
long.sim_nt <- melt(felisa_picornalike_bat_boot, id.vars = c("pointer"), measure.vars = c("OQ818339","OQ818341","HQ585111","OP880198","OQ363752","OQ363755",
                                                                                          "OQ363761","MK468720"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818339"] <- "Eidolon dupreanum picorna-like virus OQ818339"
long.sim_nt$accession[long.sim_nt$accession == "OQ818341"] <- "Eidolon dupreanum felisavirus OQ818341"
long.sim_nt$accession[long.sim_nt$accession == "HQ585111"] <- "Eptesicus fuscus picorna-like virus HQ585111"
long.sim_nt$accession[long.sim_nt$accession == "OP880198"] <- "Myotis brandtii picorna-like virus OP880198"
long.sim_nt$accession[long.sim_nt$accession == "OQ363752"] <- "Miniopterus pusillus picornavirus OQ363752"
long.sim_nt$accession[long.sim_nt$accession == "OQ363755"] <- "Miniopterus pusillus picornavirus OQ363755"
long.sim_nt$accession[long.sim_nt$accession == "OQ363761"] <- "Miniopterus pusillus picornavirus OQ363761"
long.sim_nt$accession[long.sim_nt$accession == "MK468720"] <- "Pteropus lylei picoravirales sp. MK468720"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum picorna-like virus OQ818339",
                                                                  "Eidolon dupreanum felisavirus OQ818341",
                                                                  "Eptesicus fuscus picorna-like virus HQ585111",
                                                                  "Myotis brandtii picorna-like virus OP880198",
                                                                  "Miniopterus pusillus picornavirus OQ363752",
                                                                  "Miniopterus pusillus picornavirus OQ363755",
                                                                  "Miniopterus pusillus picornavirus OQ363761",
                                                                  "Pteropus lylei picoravirales sp. MK468720"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100



## Nucleotide
title<-expression(paste("Reference: ",italic("Pteropus rufus felisavirus "), " OQ818335"))

felisapicornalike_bat_all_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,0.5,0), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,1000/3.055,2000/3.055,3000/3.055), 
  #                    labels = c(0,1000, 2000,3000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

felisapicornalike_bat_all_boot

#put gene map with PySimPlot
felispicona_all_boot<-felisapicornalike_bat_all_boot/felispicorna_all+plot_layout(nrow=2,  heights = c(2, 0.4))
felispicona_all_boot

felispicona_all_boot<-as.ggplot(felispicona_all_boot)
felispicona_all_boot



#Hepatovirus
hepatovirus_bat_all_boot <- read.csv(file = "hepatovitus_boot.csv", header = T, stringsAsFactors = F) #animo acid
head(hepatovirus_bat_all_boot)

#move to long
long.sim_nt <- melt(hepatovirus_bat_all_boot, id.vars = c("pointer"), measure.vars = c("KT452714","KX420952","MG559674","KT452730",
                                                                                       "KT452729","OM302498","KT452742"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "KT452714"] <- "Eidolon helvum hepatovirus KT452714"
long.sim_nt$accession[long.sim_nt$accession == "KX420952"] <- "Rhinopoma hardwickii hepatovirus KX420952"
long.sim_nt$accession[long.sim_nt$accession == "MG559674"] <- "Hipposideros armiger hepatovirus MG559674"
long.sim_nt$accession[long.sim_nt$accession == "KT452729"] <- "Rhinolophus landeri hepatovirus KT452729"
long.sim_nt$accession[long.sim_nt$accession == "KT452730"] <- "Coleura afra hepatovirus KT452730"
long.sim_nt$accession[long.sim_nt$accession == "OM302498"] <- "Eptesicus fuscus hepatovirus OM302498"
long.sim_nt$accession[long.sim_nt$accession == "KT452742"] <- "Miniopterus cf. manavi hepatovirus KT452742"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon helvum hepatovirus KT452714",
                                                                  "Rhinopoma hardwickii hepatovirus KX420952",
                                                                  "Hipposideros armiger hepatovirus MG559674",
                                                                  "Rhinolophus landeri hepatovirus KT452729",
                                                                  "Coleura afra hepatovirus KT452730",
                                                                  "Eptesicus fuscus hepatovirus OM302498",
                                                                  "Miniopterus cf. manavi hepatovirus KT452742"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum hepatovirus "), "OQ818337"))

hepatovirus_bat_all_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,0.5,0), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_bat_all_boot


#put gene map with PySimPlot
hep_bat_all_boot<-hepatovirus_bat_all_boot/bat_hepato_all+plot_layout(nrow=2,  heights = c(2, 0.4))
hep_bat_all_boot

hep_bat_all_boot<-as.ggplot(hep_bat_all_boot)
hep_bat_all_boot



#mischivirus
mischivirus_bat_all_boot <- read.csv(file = "mischivirus_boot.csv", header = T, stringsAsFactors = F) #Nucleotide
head(mischivirus_bat_all_boot)

#move to long
long.sim_nt <- melt(mischivirus_bat_all_boot, id.vars = c("pointer"), measure.vars = c("JQ814851","KP054273",
                                                                                       "KP054276","KP054277",
                                                                                       "KP054278",
                                                                                       "MG888045","KP100644"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "JQ814851"] <- "Miniopterus schreibersii mischivirus JQ814851"
long.sim_nt$accession[long.sim_nt$accession == "KP054273"] <- "Miniopterus schreibersii mischivirus KP054273"
long.sim_nt$accession[long.sim_nt$accession == "KP054276"] <- "Miniopterus schreibersii shanbavirus KP054276"
long.sim_nt$accession[long.sim_nt$accession == "KP054277"] <- "Miniopterus schreibersii shanbavirus KP054277"
long.sim_nt$accession[long.sim_nt$accession == "KP054278"] <- "Miniopterus schreibersii shanbavirus KP054278"
long.sim_nt$accession[long.sim_nt$accession == "MG888045"] <- "Miniopterus sp. mischivirus MG888045"
long.sim_nt$accession[long.sim_nt$accession == "KP100644"] <- "Hipposideros gigas mischivirus KP100644"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Miniopterus schreibersii mischivirus JQ814851",
                                                                  "Miniopterus schreibersii mischivirus KP054273",
                                                                  "Miniopterus schreibersii shanbavirus KP054276",
                                                                  "Miniopterus schreibersii shanbavirus KP054277",
                                                                  "Miniopterus schreibersii shanbavirus KP054278",
                                                                  "Miniopterus sp. mischivirus MG888045",
                                                                  "Hipposideros gigas mischivirus KP100644"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Pteropus rufus mischivirus "), "OQ818316"))

mischivirus_bat_all_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,0.5,0), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_bat_all_boot

#put gene map with PySimPlot
mischi_bat_all_boot<-mischivirus_bat_all_boot/bat_mischi_all+plot_layout(nrow=2,  heights = c(2, 0.4))
mischi_bat_all_boot

mischi_bat_all_boot<-as.ggplot(mischi_bat_all_boot)
mischi_bat_all_boot


#Sapelovirus
sapelovirus_bat_all_boot <- read.csv(file = "sapelovirus_boot.csv", header = T, stringsAsFactors = F) #Nucleotide
head(sapelovirus_bat_all_boot)

#move to long
long.sim_nt <- melt(sapelovirus_bat_all_boot, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","KJ641689","KJ641696",
                                                                                       "MF352431","KX644938"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_nt$accession[long.sim_nt$accession == "KJ641696"] <- "Vespertilio superans shanbavirus KJ641696"
long.sim_nt$accession[long.sim_nt$accession == "KJ641689"] <- "Myotis altarium shanbavirus KJ641689"
long.sim_nt$accession[long.sim_nt$accession == "MF352431"] <- "Myotis ricketti sapelovirus MF352431"
long.sim_nt$accession[long.sim_nt$accession == "KX644938"] <- "Eidolon helvum sapelovirus KX644938"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321",
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Vespertilio superans shanbavirus KJ641696",
                                                                  "Myotis altarium shanbavirus KJ641689",
                                                                  "Myotis ricketti sapelovirus MF352431",
                                                                  "Eidolon helvum sapelovirus KX644938"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapelovirus "), "OQ818320"))

sapelovirus_bat_all_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 7),
        legend.title = element_text(face="italic", size = 7),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,0.5,0), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  #scale_color_manual(values=colz2) + 
  scale_fill_distiller()+
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_bat_all_boot

#put gene map with PySimPlot
sapelo_bat_all_boot<-sapelovirus_bat_all_boot/bat_sapelo_all+plot_layout(nrow=2,  heights = c(2, 0.4))
sapelo_bat_all_boot

sapelo_bat_all_boot<-as.ggplot(sapelo_bat_all_boot)
sapelo_bat_all_boot

#in Fig 3
# bootscan_bat_fig<-plot_grid(felispicona_all_boot,
#                             mischi_bat_all_boot,
#                             sapelo_bat_all_boot,
#                             hep_bat_all_boot,
#                             ncol=1,
#                             labels=c("B","C","D","E"),
#                             align="hv",
#                             axis="l", label_size = 23)
# bootscan_bat_fig

bootscan_bat_fig2<-plot_grid(felispicona_all_boot,
                            mischi_bat_all_boot,
                            sapelo_bat_all_boot,
                            hep_bat_all_boot,
                            ncol=2,
                            labels=c("B","C","D","E"),
                            align="hv",
                            axis="l", label_size = 23)
bootscan_bat_fig2


#Now plot the final figure
Fig4<-plot_grid(batpicornavirales_tree,bootscan_bat_fig2,
                ncol=2,
                align="hv", axis = "b",
                rel_heights = c(1,0.1),
                rel_widths = c(1,1.3),
                label_size = 23,
                labels="A","")

Fig4

#ggsave("Fig4.pdf", width=40, height=40, units = c("in"))


