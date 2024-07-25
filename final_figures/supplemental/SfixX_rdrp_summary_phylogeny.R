#This script is for the non-collapsed and labeled version of the circular phylogeny in Fig 2, which is a master phylogeny using all the viral genera identified for the 
#RDRP gene to create an infographic

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


##Set working directory
homewd= "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/"
setwd(paste0(homewd,"/raxml_trees/polymerase_region_summary_tree"))

#load the tree and root it
tree <-  read.tree("T1.raxml.supportFBP") 
rooted.tree <- root(tree, which(tree$tip.label == "NC_001547.1"))
#take a quick look in base R
plot(rooted.tree)

#Remove root from displaying, still calculates changes correctly without it
rooted.tree<-drop.tip(rooted.tree, "NC_001547.1")

#load tree data prepared from elsewhere
dat <- read.csv(("polymerase_summary_metadata.csv"), header = T, stringsAsFactors = F)

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

dat$novel <- as.factor(dat$novel)

#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  #scale_color_manual(values=colz) + 
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

#and source only:
dat$new_label[!is.na(dat$Accession) &is.na(dat$source) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                      dat$Genus[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                      #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|", 
                                                                                      #dat$source[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                      dat$Country[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
                                                                                      dat$Collection_Date[!is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
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
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=Host), size=4,stroke=0,show.legend = T) +
  # scale_fill_manual(values=colz) +
  # scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  guides(colour = guide_legend(ncol = 1))+
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Genus="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=4, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=12), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,5))
p2

##add bootstrap values to this tree
p2.dat <- p2$data
p2.dat$Bootstrap <- NA
Bootstrap<-p2.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p2.dat$label)] <- as.numeric(p2.dat$label[(length(tree.dat$tip_label)+1):length(p2.dat$label)])#fill with label

p3 <- p2  %<+% p2.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, "cm"))
p3


# save figs
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 

ggsave(file = paste0(homewd, "/final_figures/supplemental/SfigX_rdrp_summary_phylogeny.pdf"),
       plot = p3,
       units="mm",  
       width=140, 
       height=300, 
       scale=4, 
       dpi=300)




