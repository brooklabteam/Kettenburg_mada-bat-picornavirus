rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)
library(phylotools)
library(phylobase)

homewd= "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/"

##Picornavirales all bat sequences with alignment centered around P2P3 region
setwd(paste0(homewd,"/raxml_trees/bat_picornavirales/all/bat_picorna_all_raxml/P2P3"))

#load the tree
tree <-  read.tree("T1.raxml.supportFBP") 

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("bat_picornavirales_all_P2P3_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #210
length(tree$tip.label) #210

#check subgroup names
unique(dat$Family)

colz = c("Caliciviridae" = "royalblue3",    "Picornaviridae"  = "turquoise1",
         "Coronaviridae"  = "black", "Unclassified"  = "darkorange1")

#pick order for the labels
dat$Family <- factor(dat$Family, levels = c("Picornaviridae",
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


#export collapsed 20x10 inch portrait PDF


##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Family, shape=Host), size=3,stroke=0,show.legend = T) +
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
  xlim(c(0,8))+
  geom_cladelabel(node = 324, label = "Unclassified Picornaviridae/Scotophilus kuhlii/Vietnam clade",offset=0.1, fontsize = 3, color="black") +
  geom_cladelabel(node = 225, label = "Unclassified Picornaviridae/Miniopterus/Hong Kong clade",offset=0.1, fontsize = 3, color="black") +
  geom_cladelabel(node = 258, label = "Unclassified Picornaviridae/Miniopterus pusillus/China clade",offset=0.1, fontsize = 3, color="black") +
  geom_cladelabel(node = 359, label = "Unclassified Picornaviridae/Miniopterus pusillus/China clade",offset=0.1, fontsize = 3, color="black") +
  geom_cladelabel(node = 228, label = "Parechovirus/Scotophilus kuhlii/Vietnam clade",offset=0.1, fontsize = 3, color="black") +
  geom_cladelabel(node = 293, label = "Kobuvirus/Scotophilus kuhlii/Vietnam clade",offset=0.1, fontsize = 3, color="black")
  
p2


#collapse the labeled clades
p3<-collapse(p2, 324)+geom_point2(aes(subset=(node==324)), size=3, shape=22, fill="turquoise1")
p4<-collapse(p3, 225)+geom_point2(aes(subset=(node==225)), size=3, shape=22, fill="turquoise1")
p5<-collapse(p4, 258)+geom_point2(aes(subset=(node==258)), size=3, shape=22, fill="turquoise1")
p6<-collapse(p5, 359)+geom_point2(aes(subset=(node==359)), size=3, shape=22, fill="turquoise1")
p7<-collapse(p6, 228)+geom_point2(aes(subset=(node==228)), size=3, shape=22, fill="turquoise1")
p9<-collapse(p7, 293)+geom_point2(aes(subset=(node==293)), size=3, shape=22, fill="turquoise1")
p9


##add bootstrap values to this tree
p9.dat <- p9$data
p9.dat$Bootstrap <- NA
Bootstrap<-p9.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p9.dat$label)] <- as.numeric(p9.dat$label[(length(tree.dat$tip_label)+1):length(p9.dat$label)])#fill with label

p9.1 <- p9  %<+% p9.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"))
p9.1

#export collapsed 15x15 inch PDF landscape 

