rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)
library(phylotools)
library(phylobase)
library(cowplot)
library(gggenes)
library(gridExtra)
library(grid)
library(plyr)
library(dplyr)
library(reshape2)
library(patchwork)
library(ggplotify)

###packages loaded

##Africa picornaviridae tree first
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/"

##Africa bat picornaviridae
setwd(paste0(homewd,"/raxml_trees/africa_picornaviridae/africa_picornaviridae_raxml/P2P3"))

#load the tree
tree <-  read.tree("T10.raxml.supportFBP") 

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("africa_picornaviridae_P2P3_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #356
length(tree$tip.label) #356

#check subgroup names
unique(dat$Genus)
colz = c("Cardiovirus" = "cadetblue1",    "Enterovirus"  = "cadetblue2",   "Hepatovirus"  = "cadetblue3",   
         "Kobuvirus"   = "cadetblue4" ,   "Parechovirus" = "coral1" ,"Erbovirus"  = "coral3" ,    
         "Teschovirus"  = "coral4" ,  "Sapelovirus" = "cyan1" ,   "Tremovirus"  = "cyan2" ,   
         "Anativirus"   = "cyan3" ,"Avihepatovirus"  = "cyan4","Aquamavirus"  = "darkgoldenrod1",   
         "Aphthovirus"  = "darkgoldenrod2" ,  "Senecavirus"  = "darkgoldenrod3" ,  "Cosavirus"  = "darkgoldenrod4",
         "Salivirus"   = "deepskyblue1",    "Passerivirus"  = "deepskyblue2", "Oscivirus"   = "deepskyblue3" , 
         "Unclassified picornavirus"= "deepskyblue4"  ,   "Mischivirus" = "darkorange","Pasivirus"   = "darkorange2"  ,  
         "Gallivirus"   = "darkorange3" ,  "Limnipivirus"  = "darkorange4",  "Hunnivirus"   = "firebrick1",   
         "Dicipivirus" = "firebrick2","Megrivirus"  = "firebrick3" ,   "Potamipivirus"  = "firebrick4", 
         "Sakobuvirus"  = "lightblue1" ,  "Sicinivirus"  = "lightblue2" ,  "Aalivirus"  = "lightblue3",
         "Mosavirus"  = "lightblue4" ,    "Rosavirus"  = "hotpink1" ,    "Avisivirus"   = "hotpink2",   
         "Orivirus"   = "hotpink3" ,    "Crohivirus" = "hotpink4","Torchivirus" = "indianred1" ,   
         "Bopivirus"  = "indianred3"  ,   "Malagasivirus" = "indianred4" , "Harkavirus"  = "pink1" ,   
         "Ampivirus"  = "pink2" ,"Livupivirus" = "pink3" ,   "Kunsagivirus"  = "pink4",  
         "Shanbavirus"  = "slateblue1" ,  "Rafivirus"   = "slateblue3",  "Coronavirus" ="black",  
         "Poecivirus"  = "slateblue4" ,"Rabovirus"   = "maroon1",    "Tottorivirus"  = "maroon3" , 
         "Ailurivirus" = "maroon4", "Madagascar bat kobuvirus" ="royalblue1", "Bat picornavirus"="royalblue3", "Picornavirus"="royalblue4")

#pick order for the labels
dat$Genus <- factor(dat$Genus, levels = c("Cardiovirus",    "Enterovirus",   "Hepatovirus",   
                                          "Kobuvirus",   "Parechovirus","Erbovirus",    
                                          "Teschovirus",  "Sapelovirus",   "Tremovirus",   
                                          "Anativirus","Avihepatovirus","Aquamavirus",   
                                          "Aphthovirus",  "Senecavirus",  "Cosavirus",
                                          "Salivirus",    "Passerivirus", "Oscivirus", 
                                          "Unclassified picornavirus",   "Mischivirus","Pasivirus",  
                                          "Gallivirus",  "Limnipivirus",  "Hunnivirus",   
                                          "Dicipivirus","Megrivirus",   "Potamipivirus", 
                                          "Sakobuvirus",  "Sicinivirus",  "Aalivirus",
                                          "Mosavirus",    "Rosavirus",    "Avisivirus",   
                                          "Orivirus",    "Crohivirus","Torchivirus",   
                                          "Bopivirus",   "Malagasivirus", "Harkavirus",   
                                          "Ampivirus", "Livupivirus",   "Kunsagivirus",  
                                          "Shanbavirus",  "Rafivirus", 
                                          "Poecivirus", "Rabovirus",    "Tottorivirus", 
                                          "Ailurivirus", "Madagascar bat kobuvirus", 
                                          "Bat picornavirus", "Picornavirus","Coronavirus"))   

dat$novel <- as.factor(dat$novel)


#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Genus)) +
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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Genus, Genus, novel, old_tip_label)

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
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=2, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,14))

p1

#add node shapes to represent bootstrap values
p0<-ggtree(rooted.tree)
p0.dat <- p0$data
p0.dat$Bootstrap <- NA
Bootstrap<-p0.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p0.dat$label)] <- as.numeric(p0.dat$label[(length(tree.dat$tip_label)+1):length(p0.dat$label)])#fill with label

#add bootstrap values to original plot
p1.1 <- p1  %<+% p0.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"))
p1.1


#Export uncollapsed 30x10 portrait PDF


##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=8), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,14))+
  geom_cladelabel(node = 395, label = "Apthovirus/Bos taurus clade",offset=0.1, fontsize = 3, color="black") +
  geom_cladelabel(node = 594, label = "Enterovirus/Human and non-human primate clade",offset=0.1, fontsize = 3, color="black")
p2



#collapse the labeled clades
p3<-collapse(p2, 395)+geom_point2(aes(subset=(node==395)), size=3, shape=22, fill="darkgoldenrod2")
p12<-collapse(p3, 594)+geom_point2(aes(subset=(node==594)), size=3, shape=22, fill="cadetblue2")
p12


##add bootstrap values to this tree
p12.dat <- p12$data
p12.dat$Bootstrap <- NA
Bootstrap<-p12.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p12.dat$label)] <- as.numeric(p12.dat$label[(length(tree.dat$tip_label)+1):length(p12.dat$label)])#fill with label

picornaviridae <- p12  %<+% p12.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"))
picornaviridae






##Next make the bootscan plots for part B showing recombination

##First make the genome plots
homewd="/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/genome_annotation_and_characterization/"
setwd("~/Desktop/developer/mada-bat-picornavirus/genome_annotation_and_characterization/genus_gene_maps")

#Load the gene data
ictv <- read.csv("ictv_blast_genes_trimmed.csv", header = T, stringsAsFactors = F)
ictv$gene<-factor(ictv$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                        "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                        "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Minor structural protein", "Non-structural polyprotein",
                                        "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                        "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the peptide files
ictv_pep <- read.csv("ictv_blast_peptides_trimmed.csv", header = T, stringsAsFactors = F)
ictv_pep$gene<-factor(ictv_pep$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                                "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Minor structural protein", "Non-structural polyprotein",
                                                "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))
#Load the feature file in case its needed
ictv_feat <- read.csv("ictv_blast_features_trimmed.csv", header = T, stringsAsFactors = F)
ictv_feat$gene<-factor(ictv_feat$gene, levels = c("5'UTR", "L","VP4", "VP2", "VP0", "VP3",
                                                  "VP1","VP1/2A", "2A", "2B", "2C", "3A", "3B",
                                                  "3C", "3D", "Helicase","NS4","Vpg","Pro-Pol", "Polyprotein", "Putative polyprotein", "Minor structural protein", "Non-structural polyprotein",
                                                  "Hypothetical protein", "Similar to structural polyprotein", "Structural polyprotein",
                                                  "Similar to putative polyprotein", "Similar to polyprotein","3'UTR"))

#Pick colors for genes
colz=c("5'UTR"="gold", "L"="royalblue","VP4"="paleturquoise3", "VP2"="skyblue1", "VP0"="royalblue4", "VP3"="steelblue1",
       "VP1"="cadetblue1", "VP1/2A"="cadetblue1", "2A"="orange1", "2B"="sienna2", "2C"="darkorange1", "3A"="palevioletred1", "3B"="plum",
       "3C"="rosybrown1", "3D"="pink2", "Helicase"="darkseagreen1","NS4"="darkolivegreen1","Vpg"="seagreen1","Pro-Pol"="palegreen2",
       "Polyprotein"="azure3","Putative polyprotein"="mediumorchid1", "Non-structural polyprotein"="mediumorchid4", 
       "Minor structural protein" ="slateblue3", "Structural polyprotein"="lightgoldenrod1",
       "Hypothetical protein"="darkslategrey", "Similar to structural polyprotein"="mediumpurple1", "Similar to putative polyprotein"="mediumpurple3", 
       "Similar to polyprotein"="mediumpurple4","3'UTR"="yellow")


#Plot ICTV and BLAST together plots

#Subset ICTV data by virus
ictv_batpicornavirus<-subset(ictv,molecule=="Bat picornavirus")
ictv_batpicornavirus_full<-subset(ictv_batpicornavirus,type=="full")
ictv_batpicornavirus_all<-subset(ictv_batpicornavirus,type=="all")
ictv_hepatovirus<-subset(ictv,molecule=="Hepatovirus")
ictv_kobuvirus<-subset(ictv,molecule=="Kobuvirus")
ictv_kunsagivirus<-subset(ictv,molecule=="Kunsagivirus")
ictv_mischivirus<-subset(ictv,molecule=="Mischivirus")
ictv_sapelovirus<-subset(ictv,molecule=="Sapelovirus")
ictv_sapelovirus_full<-subset(ictv_sapelovirus,type=="full")
ictv_sapelovirus_p1<-subset(ictv_sapelovirus,type=="p1")
ictv_sapelovirus_p2<-subset(ictv_sapelovirus,type=="p2")
ictv_teschovirus<-subset(ictv,molecule=="Teschovirus")

ictv_batpicornavirus_pep<-subset(ictv_pep,molecule=="Bat picornavirus")
ictv_batpicornavirus_full_pep<-subset(ictv_batpicornavirus_pep,type=="full")
ictv_batpicornavirus_all_pep<-subset(ictv_batpicornavirus_pep,type=="all")
ictv_hepatovirus_pep<-subset(ictv_pep, molecule=="Hepatovirus")
ictv_kobuvirus_pep<-subset(ictv_pep,molecule=="Kobuvirus")
ictv_kunsagivirus_pep<-subset(ictv_pep,molecule=="Kunsagivirus")
ictv_mischivirus_pep<-subset(ictv_pep,molecule=="Mischivirus")
ictv_sapelovirus_pep<-subset(ictv_pep,molecule=="Sapelovirus")
ictv_sapelovirus_pep_full<-subset(ictv_sapelovirus_pep,type=="full")
ictv_sapelovirus_pep_p1<-subset(ictv_sapelovirus_pep,type=="p1")
ictv_sapelovirus_pep_p2<-subset(ictv_sapelovirus_pep,type=="p2")
ictv_teschovirus_pep<-subset(ictv_pep,molecule=="Teschovirus")

ictv_batpicornavirus_feat<-subset(ictv_feat,molecule=="Bat picornavirus")
ictv_batpicornavirus_full_feat<-subset(ictv_batpicornavirus_feat,type=="full")
ictv_batpicornavirus_all_feat<-subset(ictv_batpicornavirus_feat,type=="all")
ictv_hepatovirus_feat<-subset(ictv_feat,molecule=="Hepatovirus")
ictv_kobuvirus_feat<-subset(ictv_feat,molecule=="Kobuvirus")
ictv_kunsagivirus_feat<-subset(ictv_feat,molecule=="Kunsagivirus")
ictv_mischivirus_feat<-subset(ictv_feat,molecule=="Mischivirus")
ictv_sapelovirus_feat<-subset(ictv_feat,molecule=="Sapelovirus")
ictv_sapelovirus_full_feat<-subset(ictv_sapelovirus_feat,type=="full")
ictv_sapelovirus_feat_p1<-subset(ictv_sapelovirus_feat,type=="p1")
ictv_sapelovirus_feat_p2<-subset(ictv_sapelovirus_feat,type=="p2")
ictv_teschovirus_feat<-subset(ictv_feat,molecule=="Teschovirus")

#plot ictv and blast plots
ictv_batpicorna_full<-ggplot(ictv_batpicornavirus_full, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_batpicornavirus_full_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                       xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_batpicornavirus_all_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(1,7700),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_batpicorna_full


ictv_hepato<-ggplot(ictv_hepatovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_hepatovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_hepatovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_hepatovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,6250),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_hepato


ictv_kun<-ggplot(ictv_kunsagivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_kunsagivirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_kunsagivirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_kunsagivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                               xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_kunsagivirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,8650),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_kun


ictv_mischi<-ggplot(ictv_mischivirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_mischivirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_mischivirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_mischivirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_mischivirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,9000),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_mischi


ictv_sapelo_p1<-ggplot(ictv_sapelovirus_p1, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_sapelovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_sapelovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_sapelovirus_pep_p1, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                                 xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_sapelovirus_feat_p1,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(810,1975),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_sapelo_p1



ictv_tescho<-ggplot(ictv_teschovirus, aes(xmin = start, xmax = end, y = molecule, fill=gene)) +
  geom_gene_arrow(arrowhead_width = grid::unit(3, "mm"),
                  arrowhead_height = grid::unit(4, "mm"),
                  arrow_body_height = grid::unit(4, "mm")) +
  # geom_feature(data=ictv_teschovirus_feat, aes(x=mid, y=molecule),
  #              feature_height = grid::unit(6,"mm"))+
  # geom_feature_label(data=ictv_teschovirus_feat, aes(x=mid, y=molecule, label=gene),
  #                    feature_height = grid::unit(6,"mm"),
  #                    label_height = grid::unit(6,"mm"))+
  geom_subgene_arrow(data = ictv_teschovirus_pep, mapping=aes(xmin = from, xmax = to, y = molecule, fill=gene,
                                                              xsubmin=from, xsubmax=to), color="black", alpha=.7,
                     arrowhead_width = grid::unit(3, "mm"),
                     arrowhead_height = grid::unit(4, "mm"),
                     arrow_body_height = grid::unit(4, "mm"))+
  geom_text(data=ictv_teschovirus_feat,mapping=aes(x = mid, y = 1.5, label = gene), size=4) +
  scale_fill_manual(values=colz)+
  theme_genes()+
  scale_x_continuous(limits=c(0,7450),expand=c(0,0))+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  xlab("Genome position") + ylab("")
ictv_tescho



##Now get the plots from bootscan

#First all bat picornavirales over 3kb
setwd("~/Desktop/developer/mada-bat-picornavirus/bootscan/output_africa")

colzpalette<-c("darkorange1","deepskyblue","goldenrod1","red3","lightpink1","slateblue1", "plum1")

#Bat picornavirus full
africa_batpicorna_full_bootscan <- read.csv(file = "africa_batpicorna_full_bootscan.csv", header = T, stringsAsFactors = F)
head(africa_batpicorna_full_bootscan)

#move to long
long.sim_nt <- melt(africa_batpicorna_full_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818325","JX437642","KF040078"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818325"] <- "Rousettus madagascariensis picornavirus OQ818325"
long.sim_nt$accession[long.sim_nt$accession == "JX437642"] <- "Homo sapiens enterovirus JX437642"
long.sim_nt$accession[long.sim_nt$accession == "KF040078"] <- "Pan troglodytes enterovirus KF040078"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rousettus madagascariensis picornavirus OQ818325", 
                                                                  "Homo sapiens enterovirus JX437642",
                                                                  "Pan troglodytes enterovirus KF040078"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Bootscan 
title<-expression(paste("Reference: ",italic("Rousettus madagascariensis picornavirus "), "OQ818328"))

batpicorna_africa_full_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  # scale_x_continuous(breaks=c(0,2000/3.055,4000/3.055,6000/3.055,8000/3.055), 
  #                    labels = c(0,2000, 4000,6000,8000),expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

batpicorna_africa_full_boot

#put gene map with PySimPlot
batpicorna_full_boot<-batpicorna_africa_full_boot/ictv_batpicorna_full+plot_layout(nrow=2,  heights = c(2, 0.30))
batpicorna_full_boot

batpicorna_full_boot<-as.ggplot(batpicorna_full_boot)
batpicorna_full_boot



#Hepatovirus
africa_hepato_bootscan <- read.csv(file = "africa_hepato_bootscan.csv", header = T, stringsAsFactors = F) #animo acid
head(africa_hepato_bootscan)

#move to long
long.sim_nt <- melt(africa_hepato_bootscan, id.vars = c("pointer"), measure.vars = c("KT452729",
                                                                                     "NC_028366","NC_038313","NC_038314"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "KT452729"] <- "Rhinolophus landeri hepatovirus KT452729"
long.sim_nt$accession[long.sim_nt$accession == "NC_028366"] <- "Eidolon helvum hepatovirus NC_028366"
long.sim_nt$accession[long.sim_nt$accession == "NC_038313"] <- "Miniopterus sp. hepatovirus NC_038313"
long.sim_nt$accession[long.sim_nt$accession == "NC_038314"] <- "Lophuromys sikapusi hepatovirus NC_038314"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Rhinolophus landeri hepatovirus KT452729",
                                                                  "Eidolon helvum hepatovirus NC_028366",
                                                                  "Miniopterus sp. hepatovirus NC_038313",
                                                                  "Lophuromys sikapusi hepatovirus NC_038314"))
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
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

hepatovirus_bat_all_boot


#put gene map with PySimPlot
hep_bat_all_boot<-hepatovirus_bat_all_boot/ictv_hepato+plot_layout(nrow=2,  heights = c(2, 0.30))
hep_bat_all_boot

hep_bat_all_boot<-as.ggplot(hep_bat_all_boot)
hep_bat_all_boot




#kunsagivirus
africa_kun_bootscan <- read.csv(file = "africa_kun_bootscan.csv", header = T, stringsAsFactors = F) #animo acid
head(africa_kun_bootscan)

#move to long
long.sim_nt <- melt(africa_kun_bootscan, id.vars = c("pointer"), measure.vars = c("NC_033818",
                                                                                  "NC_034206","HF677705"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "NC_033818"] <- "Eidolon helvum kunsagivirus NC_033818"
long.sim_nt$accession[long.sim_nt$accession == "NC_034206"] <- "Papio cynocephalus kunsagivirus NC_034206"
long.sim_nt$accession[long.sim_nt$accession == "HF677705"] <- "Hylomyscus parechovirus HF677705"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon helvum kunsagivirus NC_033818",
                                                                  "Papio cynocephalus kunsagivirus NC_034206",
                                                                  "Hylomyscus parechovirus HF677705"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum kunsagivirus "), "OQ818317"))

kunsagivirus_bat_all_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

kunsagivirus_bat_all_boot


#put gene map with PySimPlot
kun_bat_all_boot<-kunsagivirus_bat_all_boot/ictv_kun+plot_layout(nrow=2,  heights = c(2, 0.30))
kun_bat_all_boot

kun_bat_all_boot<-as.ggplot(kun_bat_all_boot)
kun_bat_all_boot



#mischivirus
africa_mischi_bootscan <- read.csv(file = "africa_mischi_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_mischi_bootscan)

#move to long
long.sim_nt <- melt(africa_mischi_bootscan, id.vars = c("pointer"), measure.vars = c("JN867757",
                                                                                     "MG888045",
                                                                                     "MN784123"))
unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "JN867757"] <- "Homo sapiens cosavirus JN867757"
long.sim_nt$accession[long.sim_nt$accession == "MG888045"] <- "Miniopterus sp. mischivirus MG888045"
long.sim_nt$accession[long.sim_nt$accession == "MN784123"] <- "Mandrillus leucophaeus cosavirus MN784123"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Homo sapiens cosavirus JN867757",
                                                                  "Miniopterus sp. mischivirus MG888045",
                                                                  "Mandrillus leucophaeus cosavirus MN784123"))
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
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

mischivirus_bat_all_boot

#put gene map with PySimPlot
mischi_bat_all_boot<-mischivirus_bat_all_boot/ictv_mischi+plot_layout(nrow=2,  heights = c(2, 0.30))
mischi_bat_all_boot

mischi_bat_all_boot<-as.ggplot(mischi_bat_all_boot)

mischi_bat_all_boot


#Sapelovirus p1
africa_sapelo_p1_bootscan <- read.csv(file = "africa_sapelo_p1_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_sapelo_p1_bootscan)

#move to long
long.sim_nt <- melt(africa_sapelo_p1_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818321","OQ818329","NC_033820",
                                                                                        "OQ818343","LC508226","OM104039"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818321"] <- "Eidolon dupreanum sapelovirus OQ818321"
long.sim_nt$accession[long.sim_nt$accession == "OQ818329"] <- "Rousettus madagascariensis sapelovirus OQ818329"
long.sim_nt$accession[long.sim_nt$accession == "NC_033820"] <- "Eidolon helvum sapelovirus NC_033820"
long.sim_nt$accession[long.sim_nt$accession == "OQ818343"] <- "Eidolon dupreanum sapelovirus OQ818343"
long.sim_nt$accession[long.sim_nt$accession == "LC508226"] <- "Porcine sapelovirus LC508226"
long.sim_nt$accession[long.sim_nt$accession == "OM104039"] <- "Porcine sapelovirus OM104039"

long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum sapelovirus OQ818321",
                                                                  "Rousettus madagascariensis sapelovirus OQ818329",
                                                                  "Eidolon helvum sapelovirus NC_033820",
                                                                  "Eidolon dupreanum sapelovirus OQ818343",
                                                                  "Porcine sapelovirus LC508226",
                                                                  "Porcine sapelovirus OM104039"))
#and plot
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Eidolon dupreanum sapelovirus "), "OQ818320"))

sapelovirus_bat_p1_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

sapelovirus_bat_p1_boot

#put gene map with PySimPlot
sapelo_bat_p1_boot<-sapelovirus_bat_p1_boot/ictv_sapelo_p1+plot_layout(nrow=2,  heights = c(2, 0.30))
sapelo_bat_p1_boot

sapelo_bat_p1_boot<-as.ggplot(sapelo_bat_p1_boot)
sapelo_bat_p1_boot


#Teschovirus
africa_tescho_bootscan <- read.csv(file = "africa_tescho_bootscan.csv", header = T, stringsAsFactors = F) #Nucleotide
head(africa_tescho_bootscan)

#move to long
long.sim_nt <- melt(africa_tescho_bootscan, id.vars = c("pointer"), measure.vars = c("OQ818318","OQ818324","OM966657","OM105029"))

unique(long.sim_nt$variable)

long.sim_nt$variable <- as.character(long.sim_nt$variable)

names(long.sim_nt)[names(long.sim_nt)=="variable"] <- "accession"

long.sim_nt$accession[long.sim_nt$accession == "OQ818318"] <- "Eidolon dupreanum teschovirus OQ818318"
long.sim_nt$accession[long.sim_nt$accession == "OQ818324"] <- "Rousettus madagascariensis teschovirus OQ818324"
long.sim_nt$accession[long.sim_nt$accession == "OM966657"] <- "Sus scrofa teschovirus OM966657"
long.sim_nt$accession[long.sim_nt$accession == "OM105029"] <- "Sus scrofa teschovirus OM105029"


long.sim_nt$accession <- factor(long.sim_nt$accession, levels = c("Eidolon dupreanum teschovirus OQ818318",
                                                                  "Rousettus madagascariensis teschovirus OQ818324",
                                                                  "Sus scrofa teschovirus OM966657",
                                                                  "Sus scrofa teschovirus OM105029"))
long.sim_nt$value[long.sim_nt$value<0] <- 0
long.sim_nt$value <- long.sim_nt$value/100

## Nucleotide
title<-expression(paste("Reference: ",italic("Rousettus madagascariensis teschovirus "), "OQ818323"))

teschovirus_bat_boot <- ggplot(long.sim_nt) + 
  geom_line(aes(x=pointer, y=value, color=accession), size=1) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(linetype = "solid", fill=NA)) + ylab("% permuted trees")+xlab("Genome position")+
  theme(panel.grid = element_blank(), strip.text = element_text(face="italic", size=12),
        strip.background = element_rect(fill="white"), 
        legend.position="top", legend.direction = "horizontal",legend.margin=margin(),
        legend.justification = "left",
        legend.text = element_text(face="italic", size = 8),
        legend.title = element_text(face="italic", size = 8),
        legend.key.height= unit(3.5, 'mm'),
        legend.key.width= unit(3.5, 'mm'),
        legend.background =element_rect(fill = alpha("white", 0)),
        axis.text = element_text(size=12), axis.title = element_text(size=12),
        plot.margin = unit(c(0,0.5,1,0.5), "cm"),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(colour = guide_legend(nrow = 3))+
  scale_color_manual(values=colzpalette) + 
  ggtitle(title)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))

teschovirus_bat_boot

#put gene map with PySimPlot
tescho_boot<-teschovirus_bat_boot/ictv_tescho+plot_layout(nrow=2,  heights = c(2, 0.30))
tescho_boot

tescho_boot<-as.ggplot(tescho_boot)
tescho_boot


##Now put the whole figure together

#bootscan fig first
bootscan<-plot_grid(mischi_bat_all_boot,
                    batpicorna_full_boot,
                    tescho_boot,
                    hep_bat_all_boot,
                    kun_bat_all_boot,
                    sapelo_bat_p1_boot,
                    ncol=2,
                    labels=c("B","C","D","E","F","G"),  label_size = 23, align = "hv", axis="b")
bootscan

#put tree side by side
fig3<-plot_grid(picornaviridae,bootscan,
                ncol=2,rel_heights = c(1,1), rel_widths = c(1,1),
                labels=c("A",""),  label_size = 23, align = "hv", axis="b")
fig3

#export 27x19inch PDF landscape


