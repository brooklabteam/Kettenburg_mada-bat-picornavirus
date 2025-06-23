#This script is for figure 2, which is individual phylogenies for each genus for which I have novel squences

rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)
library(phylotools)
library(phytools)
library(phylobase)
library(cowplot)
library(ggtreeExtra)
library(viridis)
library(ggplotify)
library(patchwork)

###packages loaded
##########################################################################################################
##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_mada-bat-picornavirus/"
setwd(paste0(homewd,"/IQtree_phylogenies/master_phylo_fig"))

#mischi
#load the tree and root it
tree <-  read.tree("mischi_align.fasta.treefile") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#compute distance between each taxa
dist.mat<-cophenetic.phylo(rooted.tree)
dist.mat<-dist.mat / max(dist.mat)
write.csv(dist.mat, "dist.mat.csv")

#load tree data prepared from elsewhere
dat <- read.csv(("mischivirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #17
length(tree$tip.label) #17

#check subgroup names
unique(dat$Species)

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Bat mischivirus 4","Bat mischivirus 5","Mischivirus A", "Mischivirus B",   "Mischivirus C",   
                                              "Mischivirus D", "Mischivirus sp.",   "Pteropus rufus mischivirus",  "Sindbis virus"))   

dat$novel <- as.factor(dat$novel)

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Species)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Species[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Species[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Species[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                                             #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])


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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Species, Species, novel, old_tip_label)

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
colz3 = c("Bat host"="#00BF7D","Non-bat host"="#00BF7D")

##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  guides(colour = guide_legend(ncol = 1))+
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  guides(colour = guide_legend(ncol = 1))+
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,7))

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
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
p1.1


##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  guides(colour = guide_legend(ncol = 1))+
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,7))
p2

#flip nodes around to show ancestral states
p2.1<-p2%>%ggtree::rotate(19)
p2.1

# p2.2<-p2.1%>%ggtree::rotate(18)
# p2.2

##add bootstrap values to this tree
p2.2dat <- p2.1$data
p2.2dat$Bootstrap <- NA
Bootstrap<-p2.2dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p2.2dat$label)] <- as.numeric(p2.2dat$label[(length(tree.dat$tip_label)+1):length(p2.2dat$label)])#fill with label

mischi <- p2.1  %<+% p2.2dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.3, "cm"))
mischi



#sapelo
#load the tree and root it
tree <-  read.tree("sapelo_align.fasta.treefile") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#compute distance between each taxa
dist.mat<-cophenetic.phylo(rooted.tree)
dist.mat<-dist.mat / max(dist.mat)
write.csv(dist.mat, "dist.mat.csv")

#load tree data prepared from elsewhere
dat <- read.csv(("sapelovirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #60
length(tree$tip.label) #60

#check subgroup names
unique(dat$Species)

#Leave this below in case you want to manually color the Species yourself
# colz = c("Bavovirus"="deepskyblue1", "Lagovirus"="darkorange",   "Minovirus"="gold1",   
#          "Nacovirus"="firebrick1",   "Nebovirus"="dodgerblue3", "Norovirus"="deeppink", "Recovirus"="aquamarine1",    
#          "Salovirus"="mediumpurple1",  "Sapovirus"="olivedrab2",   "Unclassified calicivirus"="slateblue2",   
#          "Vesivirus"="lightpink2","Alphavirus"="black")

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Bat sapelovirus", "Eidolon bat sapelovirus","Eidolon dupreanum sapelovirus 1",   "Eidolon dupreanum sapelovirus 2",   
                                              "Marmot sapelovirus 1",   "Marmot sapelovirus 2", "Pteropodidae bat sapelovirus","Rousettus bat sapelovirus","Rousettus madagascariensis sapelovirus 1", "Sapelovirus A",    
                                              "Sapelovirus B",  "Sapelovirus-like porcine picornavirus Japan", "Tasmanian devil-associated sapelovirus",  "Sindbis virus"))   

dat$novel <- as.factor(dat$novel)

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Species)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Species[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Species[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                                             #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])


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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Species, Species, novel, old_tip_label)

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
colz3 = c("Bat host"="#00BFC4","Non-bat host"="#00BFC4")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  guides(colour = guide_legend(ncol = 1))+
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  guides(colour = guide_legend(ncol = 1))+
  theme(legend.position = "none", 
        legend.direction = "horizontal",
        legend.text = element_text(size=9), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,9))

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
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        legend.key.size = unit(0.3, "cm"))
p1.1


##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  guides(colour = guide_legend(ncol = 1))+
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "horizontal",
        legend.text = element_text(size=9), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,10))+
  geom_cladelabel(node = 73, label = "Swine-hosted sapelovirus A (collapsed clade)",offset=0.1, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(60)
p2.1

p2.1<-p2.1%>%ggtree::rotate(71)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 73)+geom_point2(aes(subset=(node==73)), size=3, shape=22, fill="white")
p3

##add bootstrap values to this tree
p3.dat <- p3$data
p3.dat$Bootstrap <- NA
Bootstrap<-p3.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p3.dat$label)] <- as.numeric(p3.dat$label[(length(tree.dat$tip_label)+1):length(p3.dat$label)])#fill with label

sapelo <- p3  %<+% p3.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.3, "cm"))
sapelo



#sapo
#load the tree and root it
tree <-  read.tree("sapo_align.fasta.treefile") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#compute distance between each taxa
dist.mat<-cophenetic.phylo(rooted.tree)
dist.mat<-dist.mat / max(dist.mat)
write.csv(dist.mat, "dist.mat.csv")

#load tree data prepared from elsewhere
dat <- read.csv(("sapovirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #223
length(tree$tip.label) #223

#check subgroup names
unique(dat$Species)

#Leave this below in case you want to manually color the Species yourself
# colz = c("Bavovirus"="deepskyblue1", "Lagovirus"="darkorange",   "Minovirus"="gold1",   
#          "Nacovirus"="firebrick1",   "Nebovirus"="dodgerblue3", "Norovirus"="deeppink", "Recovirus"="aquamarine1",    
#          "Salovirus"="mediumpurple1",  "Sapovirus"="olivedrab2",   "Unclassified calicivirus"="slateblue2",   
#          "Vesivirus"="lightpink2","Alphavirus"="black")

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Bat faecal sapovirus","Bat sapovirus", "Bat sapovirus TLC34/HK", "Bat sapovirus TLC39/HK" ,"Bat sapovirus TLC58/HK",
                                              "Bat sapovirus WD3",  
                                              "Eidolon dupreanum sapovirus 1",   
                                              "Eidolon dupreanum sapovirus 2","Rousettus bat calicivirus","Rousettus madagascariensis sapovirus 2",
                                              "Rousettus madagascariensis sapovirus 4","Sapovirus sp.",  
                                              "Sapovirus Sapozj-9", "Sapovirus rat/S4-82", "Sapporo virus",    
                                              "Sindbis virus"))   

dat$novel <- as.factor(dat$novel)

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Species)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Species[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Species[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                                             #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])


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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Species, Species, novel, old_tip_label)

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
colz3 = c("Bat host"="#00B0F6","Non-bat host"="#00B0F6")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  guides(colour = guide_legend(ncol = 1))+
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  guides(colour = guide_legend(ncol = 1))+
  theme(legend.position = "bottom", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,18))

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
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        plot.margin = margin(0,0,0,0),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
p1.1


##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  guides(colour = guide_legend(ncol = 1))+
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "bottom", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,12))+
  geom_cladelabel(node = 278, label = "Human and swine-hosted Sapporo viruses (collapsed clade)",offset=0.1, fontsize=3, color="black")+
  geom_cladelabel(node = 246, label = "Hipposideros, Taphozous, and Rhinolophus bat-hosted sapoviruses (collapsed clade)",offset=0.1, fontsize=3, color="black")+
  geom_cladelabel(node = 260, label = "Myotis bat-hosted sapoviruses (collapsed clade)",offset=0.1, fontsize=3, color="black")

p2

p2.1<-p2%>%ggtree::rotate(244)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 278)+geom_point2(aes(subset=(node==278)), size=3, shape=22, fill="white")
p4.1<-collapse(p3, 260)+geom_point2(aes(subset=(node==260)), size=3, shape=22, fill="white")
p4<-collapse(p4.1, 246)+geom_point2(aes(subset=(node==246)), size=3, shape=22, fill="white")
p4


##add bootstrap values to this tree
p4.dat <- p4$data
p4.dat$Bootstrap <- NA
Bootstrap<-p4.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p4.dat$label)] <- as.numeric(p4.dat$label[(length(tree.dat$tip_label)+1):length(p4.dat$label)])#fill with label

sapo <- p4  %<+% p4.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.3, "cm"))
sapo



#kunsagi
#load the tree and root it
tree <-  read.tree("kunsagi_align.fasta.treefile") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("kunsagivirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#compute distance between each taxa
dist.mat<-cophenetic.phylo(rooted.tree)
dist.mat<-dist.mat / max(dist.mat)
write.csv(dist.mat, "dist.mat.csv")

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #5
length(tree$tip.label) #5

#check subgroup names
unique(dat$Species)

#Leave this below in case you want to manually color the Species yourself
# colz = c("Bavovirus"="deepskyblue1", "Lagovirus"="darkorange",   "Minovirus"="gold1",   
#          "Nacovirus"="firebrick1",   "Nebovirus"="dodgerblue3", "Norovirus"="deeppink", "Recovirus"="aquamarine1",    
#          "Salovirus"="mediumpurple1",  "Sapovirus"="olivedrab2",   "Unclassified calicivirus"="slateblue2",   
#          "Vesivirus"="lightpink2","Alphavirus"="black")

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Eidolon dupreanum kunsagivirus", "Kunsagivirus A",   "Kunsagivirus B",   
                                              "Kunsagivirus C", "Sindbis virus"))   

dat$novel <- as.factor(dat$novel)

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Species)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Species[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Species[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                                             #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])


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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Species, Species, novel, old_tip_label)

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
colz3 = c("Bat host"="#39B600","Non-bat host"="#39B600")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  guides(colour = guide_legend(ncol = 1))+
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  guides(colour = guide_legend(ncol = 1))+
  theme(legend.position = "none", 
        legend.direction = "horizontal",
        legend.text = element_text(size=11), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,7))

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
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
p1.1


##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  guides(colour = guide_legend(ncol = 1))+
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "horizontal",
        legend.text = element_text(size=11), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,7))
p2

##add bootstrap values to this tree
p2.dat <- p2$data
p2.dat$Bootstrap <- NA
Bootstrap<-p2.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p2.dat$label)] <- as.numeric(p2.dat$label[(length(tree.dat$tip_label)+1):length(p2.dat$label)])#fill with label

kunsagi <- p2  %<+% p2.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.3, "cm"))
kunsagi



#kobu
#load the tree and root it
tree <-  read.tree("kobu_align.fasta.treefile") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#Remove rabbit kobuvirus from displaying
rooted.tree<-drop.tip(rooted.tree, "NC_026314.1")

#compute distance between each taxa
dist.mat<-cophenetic.phylo(rooted.tree)
dist.mat<-dist.mat / max(dist.mat)
write.csv(dist.mat, "dist.mat.csv")

#load tree data prepared from elsewhere
dat <- read.csv(("kobuvirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #183
length(tree$tip.label) #183

#check subgroup names
unique(dat$Species)

#Leave this below in case you want to manually color the Species yourself
# colz = c("Bavovirus"="deepskyblue1", "Lagovirus"="darkorange",   "Minovirus"="gold1",   
#          "Nacovirus"="firebrick1",   "Nebovirus"="dodgerblue3", "Norovirus"="deeppink", "Recovirus"="aquamarine1",    
#          "Salovirus"="mediumpurple1",  "Sapovirus"="olivedrab2",   "Unclassified calicivirus"="slateblue2",   
#          "Vesivirus"="lightpink2","Alphavirus"="black")

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Aichivirus A", "Aichivirus B",   "Aichivirus C",
                                              "Aichivirus D",   "Aichivirus E", "Aichivirus F" ,"Bat kobuvirus", "Bamboo rat kobuvirus",
                                              "Canine kobuvirus", "Capreolus pygargus kobuvirus",    
                                              "Eidolon dupreanum kobuvirus",  "Kobuvirus sp.", "Marmot kobuvirus", "Mouse kobuvirus M-5/USA/2010" , 
                                              "Ovine kobuvirus",  "Sindbis virus"))   

dat$novel <- as.factor(dat$novel)

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Species)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                       dat$Species[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                       #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                       dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                       dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                       dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Species[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                      dat$Species[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                      #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                      #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                      dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                      dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                       dat$Species[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                       #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                       dat$source[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                       #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                       dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])


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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Species, Species, novel, old_tip_label)

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
colz3 = c("Bat host"="#A3A500","Non-bat host"="#A3A500")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  guides(colour = guide_legend(ncol = 1))+
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  guides(colour = guide_legend(ncol = 1))+
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,8))

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
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
p1.1


##Get the clade numbers so we can collapse unnecessary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  guides(colour = guide_legend(ncol = 1))+
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,10))+
  geom_cladelabel(node = 266, label = "Porcine and caprine-hosted Aichivirus C viruses (collapsed clade)",offset=0.1, fontsize=3, color="black") +
  geom_cladelabel(node = 330, label = "Bovine-hosted Aichivirus B viruses (collapsed clade)", offset=0.1,fontsize=3, color="black") +
  geom_cladelabel(node = 353, label = "Bovine-hosted Aichivirus C/D viruses", offset=0.05,fontsize=3, color="black") +
  geom_cladelabel(node = 194, label = "Canine, feline, and murine-hosted kobuviruses (collapsed clade)", offset=0.1,fontsize=3, color="black") +
  geom_cladelabel(node = 225, label = "Bat (Scotophilus and Rhinolophus)-hosted kobuviruses (collapsed clade)", offset=0.1,fontsize=3, color="black") +
  geom_cladelabel(node = 187, label = "Human-hosted Aichivirus A (collapsed clade)", offset=0.1,fontsize=3, color="black") +
  geom_cladelabel(node = 256, label = "Bat (Myotis and Miniopterus)-hosted Aichivirus F (collapsed clade)", offset=0.1,fontsize=3, color="black")
  
p2

p2.1<-p2%>%ggtree::rotate(182)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 266)+geom_point2(aes(subset=(node==266)), size=3, shape=22, fill="white")
p4<-collapse(p3, 330)+geom_point2(aes(subset=(node==330)), size=3, shape=22, fill="white")
p5<-collapse(p4, 353)+geom_point2(aes(subset=(node==353)), size=3, shape=22, fill="white")
p6<-collapse(p5, 187)+geom_point2(aes(subset=(node==187)), size=3, shape=22, fill="white")
p7<-collapse(p6, 225)+geom_point2(aes(subset=(node==225)), size=3, shape=22, fill="white")
p8.1<-collapse(p7, 194)+geom_point2(aes(subset=(node==194)), size=3, shape=22, fill="white")
p8<-collapse(p8.1, 256)+geom_point2(aes(subset=(node==256)), size=3, shape=22, fill="white")
p8

##add bootstrap values to this tree
p8.dat <- p8$data
p8.dat$Bootstrap <- NA
Bootstrap<-p8.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p8.dat$label)] <- as.numeric(p8.dat$label[(length(tree.dat$tip_label)+1):length(p8.dat$label)])#fill with label

kobu <- p8  %<+% p8.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.3, "cm"))
kobu


#hepato
#load the tree and root it
tree <-  read.tree("hepato_align.fasta.treefile") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#compute distance between each taxa
dist.mat<-cophenetic.phylo(rooted.tree)
dist.mat<-dist.mat / max(dist.mat)
write.csv(dist.mat, "dist.mat.csv")

#load tree data prepared from elsewhere
dat <- read.csv(("hepatovirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #129
length(tree$tip.label) #129

#check subgroup names
unique(dat$Species)

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Eidolon dupreanum hepatovirus","Eptesicus fuscus hepatovirus", 
                                              "Hepatovirus A", "Hepatovirus B",   
                                              "Hepatovirus C",   "Hepatovirus D", "Hepatovirus E", "Hepatovirus F",    
                                              "Hepatovirus G",  "Hepatovirus H",   "Hepatovirus I", "Hepatovirus sp.",
                                              "Hepatovirus sp. 'sotense'", "Hipposideros bat hepatovirus","Rodent hepatovirus", 
                                              "Taphozous bat hepatovirus" ,
                                              "Tupaia hepatovirus A",
                                              "Sindbis virus"))   

dat$novel <- as.factor(dat$novel)

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Species)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Species[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Species[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                                             #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])


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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Species, Species, novel, old_tip_label)

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
colz3 = c("Bat host"="#D89000","Non-bat host"="#D89000")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  guides(colour = guide_legend(ncol = 1))+
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  guides(colour = guide_legend(ncol = 1))+
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,8))

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
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
p1.1


##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)

#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  guides(colour = guide_legend(ncol = 1))+
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,9))+
  geom_cladelabel(node = 216, label = "Human-hosted hepatovirus A (collapsed clade)",offset=0.1, fontsize=3, color="black") +
  geom_cladelabel(node = 220, label = "Rodent,opposum, and marmot-hosted hepatovirus A, D, and F (collapsed clade)",offset=0.1, fontsize=3, color="black") +
  geom_cladelabel(node = 239, label = "Bat (Hipposideros)-hosted hepatoviruses (collapsed clade)",offset=0.1, fontsize=3, color="black") +
  geom_cladelabel(node = 233, label = "Hedgehog-hosted hepatovirus H (collapsed clade)",offset=0.1, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(129)
p2.1

p2.2<-p2.1%>%ggtree::rotate(226)
p2.2

#collapse the labeled clades
p3.1<-collapse(p2.1, 216)+geom_point2(aes(subset=(node==216)), size=4, shape=22, fill="white")
p3.2<-collapse(p3.1, 220)+geom_point2(aes(subset=(node==220)), size=4, shape=22, fill="white")
p3.3<-collapse(p3.2, 239)+geom_point2(aes(subset=(node==239)), size=4, shape=22, fill="white")
p3<-collapse(p3.3, 233)+geom_point2(aes(subset=(node==233)), size=4, shape=22, fill="white")
p3

##add bootstrap values to this tree
p3.dat <- p3$data
p3.dat$Bootstrap <- NA
Bootstrap<-p3.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p3.dat$label)] <- as.numeric(p3.dat$label[(length(tree.dat$tip_label)+1):length(p3.dat$label)])#fill with label

hepato <- p3  %<+% p3.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.3, "cm"))
hepato



#bat picorna
#load the tree and root it
tree <-  read.tree("bat_picorna_align.fasta.treefile") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#compute distance between each taxa
dist.mat<-cophenetic.phylo(rooted.tree)
dist.mat<-dist.mat / max(dist.mat)
write.csv(dist.mat, "dist.mat.csv")

#load tree data prepared from elsewhere
dat <- read.csv(("bat_picornavirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #40
length(tree$tip.label) #40

#check subgroup names
unique(dat$Species)

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Bat picornavirus BtSY4","Bat picornavirus 7","Chaerephon bat picornavirus","Shanbavirus A", "Rousettus bat picornavirus","Rousettus madagascariensis picornavirus 1", "Rousettus madagascariensis picornavirus 3",
                                              "Rousettus madagascariensis picornavirus 4","Sindbis virus"))   

dat$novel <- as.factor(dat$novel)

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Species)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Species[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Species[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                                             #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])


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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Species, Species, novel, old_tip_label)

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
colz3 = c("Bat host"="#9590FF","Non-bat host"="#9590FF")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  guides(colour = guide_legend(ncol = 1))+
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  guides(colour = guide_legend(ncol = 1))+
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,8))

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
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
p1.1


##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)

#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  guides(colour = guide_legend(ncol = 1))+
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,8))+
  geom_cladelabel(node = 46, label = "Shanbavirus A (collapsed)",offset=0.1, fontsize=3, color="black")+
  geom_cladelabel(node = 48, label = "Shanbavirus A (collapsed)",offset=0.1, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(56)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 46)+geom_point2(aes(subset=(node==46)), size=4, shape=22, fill="white")
p5<-collapse(p3, 48)+geom_point2(aes(subset=(node==48)), size=4, shape=22, fill="white")
p5

##add bootstrap values to this tree
p5.dat <- p5$data
p5.dat$Bootstrap <- NA
Bootstrap<-p5.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p5.dat$label)] <- as.numeric(p5.dat$label[(length(tree.dat$tip_label)+1):length(p5.dat$label)])#fill with label

batpicorna <- p5  %<+% p5.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.3, "cm"))
batpicorna



#tescho
#load the tree and root it
tree <-  read.tree("tescho_align.fasta.treefile") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#compute distance between each taxa
dist.mat<-cophenetic.phylo(rooted.tree)
dist.mat<-dist.mat / max(dist.mat)
write.csv(dist.mat, "dist.mat.csv")

#load tree data prepared from elsewhere
dat <- read.csv(("teschovirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #39
length(tree$tip.label) #39

#check subgroup names
unique(dat$Species)

#Leave this below in case you want to manually color the Species yourself
# colz = c("Bavovirus"="deepskyblue1", "Lagovirus"="darkorange",   "Minovirus"="gold1",   
#          "Nacovirus"="firebrick1",   "Nebovirus"="dodgerblue3", "Norovirus"="deeppink", "Recovirus"="aquamarine1",    
#          "Salovirus"="mediumpurple1",  "Sapovirus"="olivedrab2",   "Unclassified calicivirus"="slateblue2",   
#          "Vesivirus"="lightpink2","Alphavirus"="black")

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Eidolon dupreanum teschovirus 1", "Porcine teschovirus 15",   "Porcine teschovirus 16","Pteropodidae bat teschovirus",   
                                              "Rousettus madagascariensis teschovirus 1", "Rousettus bat teschovirus",  "Rousettus madagascariensis teschovirus 2", 
                                              "Teschovirus A","Teschovirus sp.","Sindbis virus"))   

dat$novel <- as.factor(dat$novel)

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Species)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Species[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Species[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                                             #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])


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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Species, Species, novel, old_tip_label)

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
colz3 = c("Bat host"="#E76BF3","Non-bat host"="#E76BF3")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  guides(colour = guide_legend(ncol = 1))+
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  guides(colour = guide_legend(ncol = 1))+
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,13))

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
  theme(legend.position = "none",
        legend.direction = "vertical",
        plot.margin = margin(0,100,0,0),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
p1.1


##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  guides(colour = guide_legend(ncol = 1))+
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,19))+
  geom_cladelabel(node = 50, label = "Porcine-hosted Teschovirus A, Porcine teschoviruses 15 and 16 (collapsed clades)",offset=0.1, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(39)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 50)+geom_point2(aes(subset=(node==50)), size=4, shape=22, fill="white")
p3

##add bootstrap values to this tree
p3.dat <- p3$data
p3.dat$Bootstrap <- NA
Bootstrap<-p3.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p3.dat$label)] <- as.numeric(p3.dat$label[(length(tree.dat$tip_label)+1):length(p3.dat$label)])#fill with label

tescho <- p3  %<+% p3.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.3, "cm"))
tescho



#cardio
#load the tree and root it
tree <-  read.tree("cardio_align.fasta.treefile") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#compute distance between each taxa
dist.mat<-cophenetic.phylo(rooted.tree)
dist.mat<-dist.mat / max(dist.mat)
write.csv(dist.mat, "dist.mat.csv")

#load tree data prepared from elsewhere
dat <- read.csv(("cardiovirus_metadata_all.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #76
length(tree$tip.label) #76

#check subgroup names
unique(dat$Species)

#Leave this below in case you want to manually color the Species yourself
# colz = c("Bavovirus"="deepskyblue1", "Lagovirus"="darkorange",   "Minovirus"="gold1",   
#          "Nacovirus"="firebrick1",   "Nebovirus"="dodgerblue3", "Norovirus"="deeppink", "Recovirus"="aquamarine1",    
#          "Salovirus"="mediumpurple1",  "Sapovirus"="olivedrab2",   "Unclassified calicivirus"="slateblue2",   
#          "Vesivirus"="lightpink2","Alphavirus"="black")

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Cardiovirus A", "Cardiovirus B",   "Cardiovirus C",   
                                              "Cardiovirus D",   "Cardiovirus E", "Cardiovirus F", "Eidolon dupreanum cardiovirus",    
                                              "Marmot cardiovirus",  "Rodent Cardiovirus",   "Sindbis virus"))   

dat$novel <- as.factor(dat$novel)

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Species)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  #scale_color_manual(values=colz) + 
  theme(legend.position = "none", legend.title = element_blank())
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Species[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Species[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Species[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
                                                                                                             #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)])


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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Species, Species, novel, old_tip_label)

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
colz3 = c("Bat host"="#F8766D","Non-bat host"="#F8766D")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  guides(colour = guide_legend(ncol = 1))+
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", family="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  guides(colour = guide_legend(ncol = 1))+
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,8))

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
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
p1.1


##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Host, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz3) +
  scale_color_manual(values=colz3)+
  scale_shape_manual(values=shapez) +
  guides(colour = guide_legend(ncol = 1))+
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05, fontface=3) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,8))+
  geom_cladelabel(node = 102, label = "Human-hosted Cardiovirus A (collapsed clades)",offset=0.1, fontsize=3, color="black") +
  geom_cladelabel(node = 96, label = "Marmot and rodent-hosted Cardiovirus E/F (collapsed clades)", offset=0.1,fontsize=3, color="black") +
  geom_cladelabel(node = 80, label = "Human and rodent-hosted Cardiovirus B/D (collapsed clades)", offset=0.1, fontsize=3, color="black")
p2

# #flip clades
# p2.1<-p2%>%ggtree::rotate(98)
# p2.1

#collapse the labeled clades
p3<-collapse(p2, 102)+geom_point2(aes(subset=(node==102)), size=4, shape=22, fill="white")
p4<-collapse(p3, 96)+geom_point2(aes(subset=(node==96)), size=4, shape=22, fill="white")
p5<-collapse(p4, 80)+geom_point2(aes(subset=(node==80)), size=4, shape=22, fill="white")
p5

##add bootstrap values to this tree
p5.dat <- p5$data
p5.dat$Bootstrap <- NA
Bootstrap<-p5.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p5.dat$label)] <- as.numeric(p5.dat$label[(length(tree.dat$tip_label)+1):length(p5.dat$label)])#fill with label

cardio <- p5  %<+% p5.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.margin = margin(0,0,0,0),
        legend.key.size = unit(0.3, "cm"))
cardio



##############################################################################################

#putting all the images together
batpicorna 
cardio #smaller
hepato 
kobu 
kunsagi #smaller
mischi
sapelo
sapo
tescho #smaller


#Plot final figure with plot_grid version 1
small_grid<-plot_grid(cardio,hepato,kobu,kunsagi,mischi, labels=c("A","B","C","D","E"),
                      rel_widths = c(1,1,1,1), rel_heights = c(0.3,0.4,0.4,0.15,0.4),
                      ncol=1, align="hv", axis="l", label_size = 23)
small_grid
small_grid<-as.ggplot(small_grid)

phylo_grid<-plot_grid(batpicorna, sapelo, tescho,sapo, labels=c("F","G","H","I"),
                      rel_widths = c(1,1,1,1), rel_heights = c(1.3,1.5,1.3,2),
                      ncol=1, align="hv", axis="l", label_size = 23)
phylo_grid
phylo_grid<-as.ggplot(phylo_grid)

final<-plot_grid(small_grid, phylo_grid, labels=c("",""),
                 rel_widths=c(1,1), rel_heights = c(1,1),
                 ncol=2,align="hv", axis="l", label_size = 23)
final

#homewd= "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/"
ggsave(file = paste0(homewd, "/final_figures/Fig2_indiv_phylogenies.pdf"),
       plot= final,
       units="mm",  
       width=250, 
       height=170, 
       scale=2, 
       dpi=500)
