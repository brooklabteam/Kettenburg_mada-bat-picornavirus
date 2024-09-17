#This script is for figure 2, which is a master phylogeny using all the viral genera identified for the 
#RDRP gene to create an infographic, then additional genus specific phylogenies to go with it. 

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
nrow(dat) #293
length(tree$tip.label) #293

#check subgroup names
unique(dat$Genus)

#pick order for the labels
dat$Genus <- factor(dat$Genus, levels = c("Cardiovirus","Hepatovirus","Kobuvirus","Kunsagivirus","Mischivirus",
                                          "Sapelovirus","Sapovirus","Teschovirus","Unclassified bat picornavirus",
                                          "Alphavirus"))   
#pick colors for virus genera
genuscolz<- c("Cardiovirus"="#F8766D","Hepatovirus"="#D89000","Kobuvirus"="#A3A500","Kunsagivirus"="#39B600","Mischivirus"="#00BF7D",
              "Sapelovirus"="#00BFC4","Sapovirus"="#00B0F6","Teschovirus"="#E76BF3","Unclassified bat picornavirus"="#9590FF",
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

novel<-data.frame(id=tree.dat$tip_label,Seq_type=tree.dat$Seq_type)
#rownames(novel) <- rooted.tree$tip.label

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
novel$Seq_type[novel$Seq_type==0] <- "Reference seq"
novel$Seq_type[novel$Seq_type==1] <- "Novel seq"
#novel$Seq_type<-as.factor(novel$Seq_type)

shapez = c("Bat host" =  17, "Non-bat host" = 19, "Reference seq"=19,"Novel seq"=17)
colz2 = c('1' =  "yellow", '0' = "white")

circ<-ggtree(rooted.tree, layout="circular")

##Get the clade numbers so we can label
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)

##base tree
p1 <- ggtree(rooted.tree, layout="fan", size=0.5) %<+% tree.dat +
  geom_tippoint(aes(color=Genus, shape=Host), size=2,stroke=0,show.legend = T) +
  #scale_color_manual(values=genuscolz)+
  scale_shape_manual(values=shapez) +
  guides(colour = "none", shape = guide_legend(order = 1))+
  theme(#legend.position = c(0.5,0.59), #keep this one in case we want the legend within the plot
        #legend.position = c(0.97,0.59), #right side
        legend.position = c(0.52,0.07),
        legend.margin = margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.25, "cm"),
        legend.direction = "vertical",
        legend.box = "horizontal") +
  xlim(c(0,15)) 
p1

library(scales)
hex_codes2 <- hue_pal()(10)      
show_col(hex_codes2)

#rotate tree a bit
p1<-rotate_tree(p1, 30)
p1

#add clade labels
p1.1 <- p1 +
  geom_cladelabel(node = 302, label = "Kobuvirus",offset=3, offset.text=0.4, fontsize=3, angle=320, hjust=0,align = TRUE, color="#A3A500") +
  geom_cladelabel(node = 472, label = "Kunsagivirus",offset=3, offset.text=0.4,  fontsize=3, angle=50,hjust=0.5,align = TRUE, color="#39B600") +
  geom_cladelabel(node = 531, label = "Cardiovirus",offset=3, offset.text=0.4,  fontsize=3, angle=30,hjust=0.5,align = TRUE, color="#F8766D") +
  geom_cladelabel(node = 556, label = "Mischivirus",offset=3, offset.text=0.4,  fontsize=3,angle=10,hjust=0.6,align = TRUE, color="#00BF7D") +
  geom_cladelabel(node = 517, label = "Teschovirus",offset=3, offset.text=0.4,  fontsize=3,angle=356,hjust=0.3,align = TRUE, color="#E76BF3") +
  geom_cladelabel(node = 329, label = "Hepatovirus",offset=3, offset.text=0.4,  fontsize=3,angle=73,hjust=0.4,align = TRUE, color="#D89000") +
  geom_cladelabel(node = 481, label = "Sapelovirus",offset=3, offset.text=0.4,  fontsize=3,angle=335,hjust=0.4,align = TRUE, color="#00BFC4") +
  geom_cladelabel(node = 500, label = "Bat picornavirus",offset=3, offset.text=0.4, angle=310,hjust=0.5, fontsize=3,align = TRUE, color="#9590FF") +
  geom_cladelabel(node = 363, label = "Sapovirus",offset=3, offset.text=0.4,angle=45, fontsize=3,align = TRUE,hjust=1, color="#00B0F6")
p1.1

##Add contig/read metadata
#attach various metadata to p1
p1.1 <- p1.1 %<+% contig

#pop contig data on top of the tree
p2<-p1.1+geom_fruit(#data=contig,
                  geom=geom_tile,
                  mapping=aes(fill=Novel_contigs),
                  pwidth=0.1,
                  offset=0.1) + 
                #scale_fill_viridis(option="G", name="Novel\ncontigs", direction = -1)
                  scale_fill_gradient(low="peachpuff1", high="orangered3")
                  
p2

#attach various metadata to p2
p2<-p2 %<+% reads

p3<-p2+new_scale_fill()+
  geom_fruit(#data=reads,
  geom=geom_tile,
  mapping=aes(fill=Novel_reads_log10),
  pwidth=0.1,
  offset=0.09)+ 
  scale_fill_gradient(low="lightblue1", high="royalblue2")+
  #scale_fill_viridis(option="B", name="Novel\nreads (log10)", direction = -1) +
 theme(
                 #legend.direction = "none",
                 legend.margin = margin(c(0,0,0,0)),
                 legend.text = element_text(size=9),
                 legend.title = element_text(size=9),
                 legend.key.size = unit(0.3, "cm"),
       plot.margin = unit(c(0, 0, 0, 0), 
                          "cm"))
p3


#attach various metadata to p3
#p3<-p3 %<+% novel

base<-p3+new_scale_fill()+
  geom_fruit(#data=reads,
    geom=geom_tile,
    mapping=aes(fill=Seq_type),
    width=0.5,
    offset=0.06)  + 
  guides(fill = guide_legend(order=2))+
  scale_fill_manual(values=c("Novel seq"="gold2","Reference seq"="grey88"), name="Seq type")+
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
colz3 = c("Bat host"="#00BF7D","Non-bat host"="#00BF7D")

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
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "horizontal",
        legend.text = element_text(size=9), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,10))+
  geom_cladelabel(node = 78, label = "Swine-hosted sapelovirus A (collapsed clade)",offset=0.1, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(60)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 78)+geom_point2(aes(subset=(node==78)), size=3, shape=22, fill="white")
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
nrow(dat) #222
length(tree$tip.label) #222

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
                                              "Eidolon dupreanum sapovirus 2","Rousettus bat calicivirus","Rousettus madagascariensis sapovirus 2","Sapovirus sp.",  
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
  theme(legend.position = "bottom",
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
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "bottom", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,12))+
  geom_cladelabel(node = 276, label = "Human and swine-hosted Sapporo viruses (collapsed clade)",offset=0.1, fontsize=3, color="black")+
  geom_cladelabel(node = 245, label = "Myotis bat-hosted sapoviruses (collapsed clade)",offset=0.1, fontsize=3, color="black")

p2

p2.1<-p2%>%ggtree::rotate(243)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 276)+geom_point2(aes(subset=(node==276)), size=3, shape=22, fill="white")
p4<-collapse(p3, 245)+geom_point2(aes(subset=(node==245)), size=3, shape=22, fill="white")
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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,10))+
  geom_cladelabel(node = 296, label = "Porcine and caprine-hosted Aichivirus C viruses (collapsed clade)",offset=0.1, fontsize=3, color="black") +
  geom_cladelabel(node = 273, label = "Bovine-hosted Aichivirus B viruses (collapsed clade)", offset=0.1,fontsize=3, color="black") +
  #geom_cladelabel(node = 241, label = "Aichivirus A (collapsed)", offset=0.05,fontsize=3, color="black") +
  geom_cladelabel(node = 255, label = "Human-hosted Aichivirus A (collapsed clade)", offset=0.1,fontsize=3, color="black") +
  geom_cladelabel(node = 198, label = "Bat (Scotophilus and Rhinolophus)-hosted kobuvirused (collapsed clade)", offset=0.1,fontsize=3, color="black") +
  geom_cladelabel(node = 225, label = "Canine, feline, and murine-hosted kobuviruses (collapsed clade)", offset=0.1,fontsize=3, color="black") +
  geom_cladelabel(node = 185, label = "Bat (Myotis and Miniopterus)-hosted Aichivirus F (collapsed clade)", offset=0.1,fontsize=3, color="black") +
  geom_cladelabel(node = 264, label = "Bovine-hosted Aichivirus D viruses (collapsed clade)", offset=0.1,fontsize=3, color="black")
  
p2

p2.1<-p2%>%ggtree::rotate(193)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 296)+geom_point2(aes(subset=(node==296)), size=3, shape=22, fill="white")
p4<-collapse(p3, 273)+geom_point2(aes(subset=(node==273)), size=3, shape=22, fill="white")
#p5<-collapse(p4, 241)+geom_point2(aes(subset=(node==241)), size=3, shape=22, fill="white")
p6<-collapse(p4, 255)+geom_point2(aes(subset=(node==255)), size=3, shape=22, fill="white")
p7<-collapse(p6, 198)+geom_point2(aes(subset=(node==198)), size=3, shape=22, fill="white")
p8.1<-collapse(p7, 225)+geom_point2(aes(subset=(node==225)), size=3, shape=22, fill="white")
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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,9))+
  geom_cladelabel(node = 185, label = "Human-hosted hepatovirus A (collapsed clade)",offset=0.1, fontsize=3, color="black") +
  geom_cladelabel(node = 218, label = "Rodent,opposum, and marmot-hosted hepatovirus A, D, and F (collapsed clade)",offset=0.1, fontsize=3, color="black") +
  geom_cladelabel(node = 204, label = "Bat (Hipposideros)-hosted hepatoviruses (collapsed clade)",offset=0.1, fontsize=3, color="black") +
  geom_cladelabel(node = 197, label = "Hedgehog-hosted hepatovirus H (collapsed clade)",offset=0.1, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(129)
p2.1

p2.2<-p2.1%>%ggtree::rotate(190)
p2.2

p2.3<-p2.2%>%ggtree::rotate(194)
p2.3

#collapse the labeled clades
p3.1<-collapse(p2.3, 185)+geom_point2(aes(subset=(node==185)), size=4, shape=22, fill="white")
p3.2<-collapse(p3.1, 218)+geom_point2(aes(subset=(node==218)), size=4, shape=22, fill="white")
p3.3<-collapse(p3.2, 204)+geom_point2(aes(subset=(node==204)), size=4, shape=22, fill="white")
p3<-collapse(p3.3, 197)+geom_point2(aes(subset=(node==197)), size=4, shape=22, fill="white")
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
nrow(dat) #39
length(tree$tip.label) #39

#check subgroup names
unique(dat$Species)

#pick order for the labels
dat$Species <- factor(dat$Species, levels = c("Bat picornavirus BtSY4","Bat picornavirus 7","Chaerephon bat picornavirus","Shanbavirus A", "Rousettus bat picornavirus","Rousettus madagascariensis picornavirus 1", "Rousettus madagascariensis picornavirus 3",
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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,8))+
  geom_cladelabel(node = 45, label = "Shanbavirus A (collapsed)",offset=0.1, fontsize=3, color="black")+
  geom_cladelabel(node = 48, label = "Shanbavirus A (collapsed)",offset=0.1, fontsize=3, color="black")+
  geom_cladelabel(node = 70, label = "Shanbavirus A (collapsed)",offset=0.1, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(58)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 45)+geom_point2(aes(subset=(node==45)), size=4, shape=22, fill="white")
p4<-collapse(p3, 48)+geom_point2(aes(subset=(node==48)), size=4, shape=22, fill="white")
p5<-collapse(p4, 70)+geom_point2(aes(subset=(node==70)), size=4, shape=22, fill="white")
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
nrow(dat) #38
length(tree$tip.label) #38

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
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,10))+
  geom_cladelabel(node = 51, label = "Porcine-hosted Teschovirus A, Porcine teschoviruses 15 and 16 (collapsed clades)",offset=0.1, fontsize=3, color="black")
p2

p2.1<-p2%>%ggtree::rotate(38)
p2.1

#collapse the labeled clades
p3<-collapse(p2.1, 51)+geom_point2(aes(subset=(node==51)), size=4, shape=22, fill="white")
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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Species="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-2, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,8))+
  geom_cladelabel(node = 95, label = "Human-hosted Cardiovirus A (collapsed clades)",offset=0.1, fontsize=3, color="black") +
  geom_cladelabel(node = 107, label = "Marmot and rodent-hosted Cardiovirus E/F (collapsed clades)", offset=0.1,fontsize=3, color="black") +
  geom_cladelabel(node = 109, label = "Human and rodent-hosted Cardiovirus B/D (collapsed clades)", offset=0.1, fontsize=3, color="black")
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


#Plot final figure with plot_grid version 1
phylo_grid<-plot_grid( mischi, batpicorna, sapelo, tescho,sapo, labels=c("F","G","H","I","J"),
                      rel_widths = c(1,1,1,1,1), rel_heights = c(0.9,1,1,0.6,1.7),
                      ncol=1, align="hv", axis="l", label_size = 23)
phylo_grid
phylo_grid<-as.ggplot(phylo_grid)

small_grid<-plot_grid(cardio,hepato,kobu,kunsagi, labels=c("B","C","D","E"),
                      rel_widths = c(1,1,1), rel_heights = c(0.6,1,1,0.3),
                      ncol=1, align="hv", axis="l", label_size = 23)
small_grid
small_grid<-as.ggplot(small_grid)

leftside<-plot_grid(base, small_grid, labels=c("A",""),
               rel_widths=c(2,1), rel_heights=c(1,1),
               ncol=1, align="h", axis="l", label_size=23)
leftside
leftside<-as.ggplot(leftside)


final<-plot_grid(leftside, phylo_grid, labels=c("",""),
                 rel_widths=c(1,1), rel_heights = c(1,1),
                 ncol=2,align="hv", axis="l", label_size = 23)
final



#homewd= "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/"
ggsave(file = paste0(homewd, "/final_figures/Fig2_summary_phylogeny.pdf"),
       plot= final,
       units="mm",  
       width=235, 
       height=230, 
       scale=2, 
       dpi=500)

