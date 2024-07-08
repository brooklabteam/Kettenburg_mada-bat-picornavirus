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

###packages loaded


#all tree code
##########################################################################################################
##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/"
setwd(paste0(homewd,"/raxml_trees/master_phylo_fig"))

#load the tree and root it
tree <-  read.tree("circ.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#load tree data prepared from elsewhere
dat <- read.csv(("polymerase_summary_metadata_circ.csv"), header = T, stringsAsFactors = F)

head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #287
length(tree$tip.label) #287

#check subgroup names
unique(dat$Genus)

#pick order for the labels
dat$Genus <- factor(dat$Genus, levels = c("Cardiovirus","Hepatovirus","Kobuvirus","Kunsagivirus","Mischivirus",
                                          "Sapelovirus","Sapovirus","Shanbavirus","Teschovirus","Unclassified picornavirus",
                                          "Alphavirus"))   
#pick colors for virus genera
genuscolz<- c("Cardiovirus"="#0A9F9D","Hepatovirus"="#CEB175","Kobuvirus"="#E54E21","Kunsagivirus"="#6C8645","Mischivirus"="#C18748",
              "Sapelovirus"="#C52E19","Sapovirus"="#AF4E24","Shanbavirus"="#54D8B1","Teschovirus"="#b67c3b","Unclassified picornavirus"="#175149",
              "Alphavirus"="black")

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  scale_color_manual(values=genuscolz) + 
  theme(legend.position = "none", legend.title = element_blank())
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
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Genus[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                             dat$source[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Country[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                             dat$Collection_Date[!is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|", 
#                                                                                                             dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                             dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) &!is.na(dat$Country)])
# 
#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            dat$Genus[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                                            #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                                            dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[ !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
                                                                                                             dat$Genus[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
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

#making separate datasets for heatmaps to go next to tree
extra<-data.frame(id=tree.dat$tip_label,Novel_contigs=tree.dat$num_genus_novel_contigs,
                  Novel_reads_log10=tree.dat$num_genus_novel_reads_log,Region=tree.dat$region,
                  Host_class=tree.dat$host_class)
#rownames(extra) <- rooted.tree$tip.label

contig<-data.frame(id=tree.dat$tip_label,Novel_contigs=tree.dat$num_genus_novel_contigs)
#rownames(contig) <- rooted.tree$tip.label

reads<-data.frame(id=tree.dat$tip_label,Novel_reads_log10=tree.dat$num_genus_novel_reads_log)
#rownames(reads) <- rooted.tree$tip.label

region<-data.frame(id=tree.dat$tip_label,Region=tree.dat$region)
#rownames(region) <- rooted.tree$tip.label

host<-data.frame(id=tree.dat$tip_label,Host_class=tree.dat$host_class)
#rownames(host) <- rooted.tree$tip.label

#make real tree.dat file for putting the tip points on the tree
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, source, Country, Collection_Date, Genus, Seq_type, old_tip_label)

rooted.tree$tip.label <- tree.dat$tip_label

head(tree.dat)
#verify that the old labels match the new

cbind(tree.dat$old_tip_label, rooted.tree$tip.label)

#check out the labels
tree.dat$tip_label#all look good

#assign some stuff for shapes
tree.dat$Host[tree.dat$Host==0] <- "Non-bat host"
tree.dat$Host[tree.dat$Host==1] <- "Bat host"
tree.dat$Host <- as.factor(tree.dat$Host)
tree.dat$Seq_type[tree.dat$Seq_type==0] <- "Reference seq"
tree.dat$Seq_type[tree.dat$Seq_type==1] <- "Novel seq"
tree.dat$Seq_type<-as.factor(tree.dat$Seq_type)

shapez = c("Bat host" =  17, "Non-bat host" = 19, "Reference seq"=19,"Novel seq"=17)
colz2 = c('1' =  "yellow", '0' = "white")

circ<-ggtree(rooted.tree, layout="circular")

##Get the clade numbers so we can label
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)

##base tree
p1 <- ggtree(rooted.tree, layout = "circular") %<+% tree.dat +
  geom_tippoint(aes(color=Genus, shape=Seq_type), size=2,stroke=0,show.legend = T) +
  scale_color_manual(values=genuscolz)+
  scale_shape_manual(values=shapez) +
  guides(colour = "none")+
  theme(legend.position = c(0.305,0.505),
        legend.direction = "vertical",
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        #legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.25, "cm")) +
  xlim(c(0,5)) 
p1

#rotate tree a bit
p1<-rotate_tree(p1, 90)
p1

#add clade labels
p1.1 <- p1 +
  geom_cladelabel(node = 290, label = "Kobuvirus",offset.text=0.1, fontsize=3, angle=315, hjust=0,align = TRUE, color="#E54E21") +
  geom_cladelabel(node = 288, label = "Kunsagivirus",offset.text=0.1, fontsize=3, angle=280,hjust=0.5,align = TRUE, color="#6C8645") +
  geom_cladelabel(node = 352, label = "Cardiovirus",offset.text=0.1, fontsize=3, angle=83,hjust=0.6,align = TRUE, color="#0A9F9D") +
  geom_cladelabel(node = 377, label = "Mischivirus",offset.text=0.1, fontsize=3,angle=65,hjust=0.4,align = TRUE, color="#C18748") +
  geom_cladelabel(node = 382, label = "Teschovirus",offset.text=0.1, fontsize=3,angle=51,hjust=0.7,align = TRUE, color="#b67c3b") +
  geom_cladelabel(node = 535, label = "Hepatovirus",offset.text=0.1, fontsize=3,angle=25,hjust=0.4,align = TRUE, color="#CEB175") +
  geom_cladelabel(node = 530, label = "Shanbavirus",offset.text=0.1, fontsize=3,angle=358,hjust=0.4,align = TRUE, color="#54D8B1") +
  geom_cladelabel(node = 503, label = "Sapelovirus",offset.text=0.1, fontsize=3,angle=340,hjust=0.4,align = TRUE, color="#C52E19") +
  geom_cladelabel(node = 501, label = "Bat picornavirus",offset.text=0.1,angle=320,hjust=0.5, fontsize=3,align = TRUE, color="#175149") +
  geom_cladelabel(node = 396, label = "Sapovirus",offset.text=0.1,angle=75, fontsize=3,align = TRUE,hjust=1, color="#AF4E24")
p1.1


##Add contig/read metadata
#attach various metadata to p1
p1.1 <- p1.1 %<+% contig

#pop contig data on top of the tree
p2<-p1.1+geom_fruit(#data=contig,
                  geom=geom_tile,
                  mapping=aes(fill=Novel_contigs),
                  width=0.3,
                  offset=0.15) + 
                scale_fill_viridis(option="G")
                  
p2

#attach various metadata to p2
p2<-p2 %<+% reads

base<-p2+new_scale_fill()+
  geom_fruit(#data=reads,
  geom=geom_tile,
  mapping=aes(fill=Novel_reads_log10),
  width=0.3,
  offset=0.08
)  + scale_fill_viridis(option="B") +
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))
 theme(
                 #legend.direction = "none",
                 legend.text = element_text(size=9),
                 legend.title = element_text(size=9),
                 legend.key.size = unit(0.3, "cm"),
       plot.margin = unit(c(0, 0, 0, 0), 
                          "cm"))
base



#mischi
#load the tree and root it
tree <-  read.tree("mischi.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

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
colz3 = c("Bat host"="black","Non-bat host"="grey50")

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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,7))
p2

#flip nodes around to show ancestral states
p2.1<-p2%>%ggtree::rotate(20)
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
        legend.key.size = unit(0.3, "cm"))
mischi



#sapelo
#load the tree and root it
tree <-  read.tree("sapelo.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#load tree data prepared from elsewhere
dat <- read.csv(("sapelovirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #57
length(tree$tip.label) #57

#check subgroup names
unique(dat$Species)

#Leave this below in case you want to manually color the Species yourself
# colz = c("Bavovirus"="deepskyblue1", "Lagovirus"="darkorange",   "Minovirus"="gold1",   
#          "Nacovirus"="firebrick1",   "Nebovirus"="dodgerblue3", "Norovirus"="deeppink", "Recovirus"="aquamarine1",    
#          "Salovirus"="mediumpurple1",  "Sapovirus"="olivedrab2",   "Unclassified calicivirus"="slateblue2",   
#          "Vesivirus"="lightpink2","Alphavirus"="black")

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Bat sapelovirus", "Eidolon dupreanum sapelovirus 1",   "Eidolon dupreanum sapelovirus 2",   
                                              "Marmot sapelovirus 1",   "Marmot sapelovirus 2", "Pteropodidae bat sapelovirus","Rousettus madagascariensis sapelovirus 1", "Sapelovirus A",    
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
colz3 = c("Bat host"="black","Non-bat host"="grey50")


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
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "horizontal",
        legend.text = element_text(size=9), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,10))+
  geom_cladelabel(node = 70, label = "Sapelovirus A/B (collapsed)",offset=0.05, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(57)
p2.1

p2.2<-p2.1%>%ggtree::rotate(65)
p2.2

p2.3<-p2.2%>%ggtree::rotate(58)
p2.3

#collapse the labeled clades
p3<-collapse(p2.3, 70)+geom_point2(aes(subset=(node==70)), size=3, shape=22, fill="white")
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
        legend.key.size = unit(0.3, "cm"))
sapelo



#sapo
#load the tree and root it
tree <-  read.tree("sapo.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#load tree data prepared from elsewhere
dat <- read.csv(("sapovirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #218
length(tree$tip.label) #218

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
                                              "Eidolon dupreanum sapovirus 2","Rousettus madagascariensis sapovirus 2","Sapovirus sp.",  
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
colz3 = c("Bat host"="black","Non-bat host"="grey50")


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
  theme(legend.position = "right", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,9))+
  geom_cladelabel(node = 336, label = "Sapporo virus (collapsed)",offset=0.05, fontsize=3, color="black") +
  geom_cladelabel(node = 277, label = "Sapporo virus/Bat sapovirus (collapsed)", offset=0.05,fontsize=3, color="black") +
  geom_cladelabel(node = 315, label = "Sapporo virus (collapsed)", offset=0.05, fontsize=3, color="black")+
  geom_cladelabel(node = 282, label = "Bat sapovirus (collapsed)", offset=0.05, fontsize=3, color="black")

p2

p2.1<-p2%>%ggtree::rotate(279)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 336)+geom_point2(aes(subset=(node==336)), size=3, shape=22, fill="white")
p4<-collapse(p3, 277)+geom_point2(aes(subset=(node==277)), size=3, shape=22, fill="white")
p6.1<-collapse(p4, 315)+geom_point2(aes(subset=(node==315)), size=3, shape=22, fill="white")
p6<-collapse(p6.1, 282)+geom_point2(aes(subset=(node==282)), size=3, shape=22, fill="white")
p6


##add bootstrap values to this tree
p6.dat <- p6$data
p6.dat$Bootstrap <- NA
Bootstrap<-p6.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p6.dat$label)] <- as.numeric(p6.dat$label[(length(tree.dat$tip_label)+1):length(p6.dat$label)])#fill with label

sapo <- p6  %<+% p6.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
sapo



#kunsagi
#load the tree and root it
tree <-  read.tree("kunsagi.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("kunsagivirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

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
colz3 = c("Bat host"="black","Non-bat host"="grey50")


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
  geom_treescale(fontsize=4, x=0,y=0, linesize = .5) +
  guides(colour = guide_legend(ncol = 1))+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
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
  theme(legend.position = "bottom",
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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=0, linesize = .5) +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.text = element_text(size=12), 
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
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        legend.key.size = unit(0.3, "cm"))
kunsagi



#kobu
#load the tree and root it
tree <-  read.tree("kobu.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#Remove rabbit kobuvirus from displaying
rooted.tree<-drop.tip(rooted.tree, "NC_026314.1")

#load tree data prepared from elsewhere
dat <- read.csv(("kobuvirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #183
length(tree$tip.label) #182

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
colz3 = c("Bat host"="black","Non-bat host"="grey50")


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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,10))+
  geom_cladelabel(node = 296, label = "Aichivirus C (collapsed)",offset=0.05, fontsize=3, color="black") +
  geom_cladelabel(node = 273, label = "Aichivirus B (collapsed)", offset=0.05,fontsize=3, color="black") +
  geom_cladelabel(node = 241, label = "Aichivirus A (collapsed)", offset=0.05,fontsize=3, color="black") +
  geom_cladelabel(node = 255, label = "Aichivirus A (collapsed)", offset=0.05,fontsize=3, color="black") +
  geom_cladelabel(node = 201, label = "Scotophilus kobuvirus (collapsed)", offset=0.05,fontsize=3, color="black") +
  geom_cladelabel(node = 227, label = "Canine kobuvirus (collapsed)", offset=0.05,fontsize=3, color="black") +
  geom_cladelabel(node = 185, label = "Aichivirus F (collapsed)", offset=0.05,fontsize=3, color="black") +
  geom_cladelabel(node = 264, label = "Aichivirus D (collapsed)", offset=0.05,fontsize=3, color="black")
  
p2

p2.1<-p2%>%ggtree::rotate(193)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 296)+geom_point2(aes(subset=(node==296)), size=3, shape=22, fill="white")
p4<-collapse(p3, 273)+geom_point2(aes(subset=(node==273)), size=3, shape=22, fill="white")
p5<-collapse(p4, 241)+geom_point2(aes(subset=(node==241)), size=3, shape=22, fill="white")
p6<-collapse(p5, 255)+geom_point2(aes(subset=(node==255)), size=3, shape=22, fill="white")
p7<-collapse(p6, 201)+geom_point2(aes(subset=(node==201)), size=3, shape=22, fill="white")
p8.1<-collapse(p7, 227)+geom_point2(aes(subset=(node==227)), size=3, shape=22, fill="white")
p8.2<-collapse(p8.1, 185)+geom_point2(aes(subset=(node==185)), size=3, shape=22, fill="white")
p8<-collapse(p8.2, 264)+geom_point2(aes(subset=(node==264)), size=3, shape=22, fill="white")
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
        legend.key.size = unit(0.3, "cm"))
kobu


#hepato
#load the tree and root it
tree <-  read.tree("hepato.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

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
colz3 = c("Bat host"="black","Non-bat host"="grey50")


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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,9))+
  geom_cladelabel(node = 189, label = "Hepatovirus A/D/F/I (collapsed)",offset=0.05, fontsize=3, color="black") +
  geom_cladelabel(node = 213, label = "Bat hepatovirus (collapsed)",offset=0.05, fontsize=3, color="black") +
  geom_cladelabel(node = 203, label = "Bat hepatovirus (collapsed)",offset=0.05, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(129)
p2.1

p2.2<-p2.1%>%ggtree::rotate(190)
p2.2

p2.3<-p2.2%>%ggtree::rotate(194)
p2.3

#collapse the labeled clades
p3.1<-collapse(p2.3, 189)+geom_point2(aes(subset=(node==189)), size=4, shape=22, fill="white")
p3.2<-collapse(p3.1, 213)+geom_point2(aes(subset=(node==213)), size=4, shape=22, fill="white")
p3<-collapse(p3.2, 203)+geom_point2(aes(subset=(node==203)), size=4, shape=22, fill="white")
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
        legend.key.size = unit(0.3, "cm"))
hepato



#bat picorna
#load the tree and root it
tree <-  read.tree("batpicorna.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#load tree data prepared from elsewhere
dat <- read.csv(("bat_picornavirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #31
length(tree$tip.label) #31

#check subgroup names
unique(dat$Species)

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Shanbavirus A", "Rousettus madagascariensis picornavirus 1", "Rousettus madagascariensis picornavirus 3",
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
colz3 = c("Bat host"="black","Non-bat host"="grey50")


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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,8))+
  geom_cladelabel(node = 38, label = "Shanbavirus A (collapsed)",offset=0.05, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(31)
p2.1

p2.2<-p2.1%>%ggtree::rotate(39)
p2.2

#collapse the labeled clades
p3<-collapse(p2.2, 38)+geom_point2(aes(subset=(node==38)), size=4, shape=22, fill="white")
p3

##add bootstrap values to this tree
p3.dat <- p3$data
p3.dat$Bootstrap <- NA
Bootstrap<-p3.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p3.dat$label)] <- as.numeric(p3.dat$label[(length(tree.dat$tip_label)+1):length(p3.dat$label)])#fill with label

batpicorna <- p3  %<+% p3.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
batpicorna



#tescho
#load the tree and root it
tree <-  read.tree("tescho.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#load tree data prepared from elsewhere
dat <- read.csv(("teschovirus_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #36
length(tree$tip.label) #36

#check subgroup names
unique(dat$Species)

#Leave this below in case you want to manually color the Species yourself
# colz = c("Bavovirus"="deepskyblue1", "Lagovirus"="darkorange",   "Minovirus"="gold1",   
#          "Nacovirus"="firebrick1",   "Nebovirus"="dodgerblue3", "Norovirus"="deeppink", "Recovirus"="aquamarine1",    
#          "Salovirus"="mediumpurple1",  "Sapovirus"="olivedrab2",   "Unclassified calicivirus"="slateblue2",   
#          "Vesivirus"="lightpink2","Alphavirus"="black")

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Eidolon dupreanum teschovirus 1", "Porcine teschovirus 15",   "Porcine teschovirus 16","Pteropodidae bat teschovirus",   
                                              "Rousettus madagascariensis teschovirus 1",   "Rousettus madagascariensis teschovirus 2", 
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
colz3 = c("Bat host"="black","Non-bat host"="grey50")


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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,10))+
  geom_cladelabel(node = 44, label = "Teschovirus A (collapsed)",offset=0.05, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(36)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 44)+geom_point2(aes(subset=(node==44)), size=4, shape=22, fill="white")
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
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
tescho



#cardio
#load the tree and root it
tree <-  read.tree("cardio.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

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
colz3 = c("Bat host"="black","Non-bat host"="grey50")


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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-1, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,8))+
  geom_cladelabel(node = 95, label = "Cardiovirus A (collapsed)",offset=0.05, fontsize=3, color="black") +
  geom_cladelabel(node = 107, label = "Cardiovirus E/F (collapsed)", offset=0.05,fontsize=3, color="black") +
  geom_cladelabel(node = 109, label = "Cardiovirus B/D/Rodent cardiovirus (collapsed)", offset=0.05, fontsize=3, color="black")
p2

#flip clades
p2.1<-p2%>%ggtree::rotate(98)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 95)+geom_point2(aes(subset=(node==95)), size=4, shape=22, fill="white")
p4<-collapse(p3, 107)+geom_point2(aes(subset=(node==107)), size=4, shape=22, fill="white")
p5<-collapse(p4, 109)+geom_point2(aes(subset=(node==109)), size=4, shape=22, fill="white")
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
        legend.key.size = unit(0.3, "cm"))
cardio



##############################################################################################

#putting all the images together
base #base figure
batpicorna 
cardio #smaller
hepato 
kobu 
kunsagi #smaller
mischi
sapelo
sapo
tescho #smaller


smallphylo <-plot_grid(sapelo, tescho, cardio, labels=c("X","X","X"), 
                       rel_widths = c(1,1,1), rel_heights = c(1,1,1),
                       ncol=1, align="hv", axis="l", label_size = 23)
smallphylo
smallphylo<-as.ggplot(smallphylo)


top_grid<-plot_grid(base,smallphylo, labels=c("X",""),
                    rel_widths = c(3,1.5), rel_heights = c(1,1),
                    ncol=2, align="hv", axis="l", label_size = 23)
top_grid
top_grid<-as.ggplot(top_grid)


doublegrid<-plot_grid(batpicorna,hepato,kobu,mischi, sapo, kunsagi, labels=c("X","X","X","X","X","X"),
                      rel_widths = c(1,1,1,1,1,1), rel_heights = c(1,1,1,1,1,2),
                      nrow = 3, align="hv", axis="l", label_size = 23)
doublegrid
doublegrid<-as.ggplot(doublegrid)


final<- plot_grid(top_grid,doublegrid, labels=c("",""),
                  rel_widths = c(1,1), rel_heights = c(1,1),
                  ncol=1, align="hv", axis="b", label_size = 23)
final







