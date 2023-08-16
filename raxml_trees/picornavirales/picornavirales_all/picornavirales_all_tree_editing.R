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

##Picornavirales all full and partial sequences >3kb
setwd(paste0(homewd,"/raxml_trees/picornavirales/picornavirales_all"))

#load the tree
tree <-  read.tree("T10.raxml.supportFBP") 

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("picornavirales_all_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #416
length(tree$tip.label) #416

#check subgroup names
unique(dat$Family)

colz = c("Caliciviridae" = "royalblue3",    "Picornaviridae"  = "turquoise1",   "Secoviridae"  = "goldenrod1",   "Dicistroviridae"   = "dodgerblue1" ,   "Iflaviridae" = "firebrick1" ,
         "Polycipiviridae"  = "lightpink1" ,    "Marnaviridae"  = "hotpink1" ,  "Solinviviridae" = "lightskyblue" ,   "Coronaviridae"  = "black", "Unclassified"  = "darkorange1")

#pick order for the labels
dat$Family <- factor(dat$Family, levels = c("Dicistroviridae",    "Secoviridae",   "Marnaviridae",   "Picornaviridae",   "Solinviviridae",
                                            "Caliciviridae",    "Iflaviridae", "Polycipiviridae",  "Unclassified","Coronaviridae"))   

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

#and accession number only
#check:
dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) &!is.na(dat$source) & !is.na(dat$Country)] #none with this condition



#now NA in Isolate + source
#first check if there are any
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)] #no
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)]<- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                        dat$Family[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                        #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
#                                                                                                        #dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|",
#                                                                                                        dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                        dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])
#Isolate + accession only
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) &!is.na(dat$Country)] #no

#Isolate + Country only
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)] #no
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|", 
#                                                                                                                 dat$Family[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|", 
#                                                                                                           #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
#                                                                                                           dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|",
#                                                                                                           #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                           dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)])

#now look at the other 2-way NAs: accession + source, accession + Country, source+Country
dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$source) & !is.na(dat$Country)] #none

dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) & is.na(dat$Country)] #none
dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)] #none

#replace
# dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)]<- paste(dat$Accession[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|", 
#                                                                                                                 dat$Family[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|", 
#                                                                                                                 dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|", 
#                                                                                                                 #dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|",
#                                                                                                                 #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                                 dat$Collection_Date[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)])

#and 3-way NAs: Isolate + accession + source
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$source) & !is.na(dat$Country)] #none
#Isolate + accession + Country
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) & is.na(dat$Country)] #none
#Isolate source Country
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)] #none
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|",
#                                                                                                                 dat$Family[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|",
#                                                                                                                 #dat$Isolate[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|",
#                                                                                                                 #dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|",
#                                                                                                                 #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) &!is.na(dat$Country)], "|",
#                                                                                                                 dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)])
# #accession Country source
dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)] #none
#and all four
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)] #none

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
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"))
p1.1

# p1.2<-p1.1%>%ggtree::rotate(node=570)
# p1.2

##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Family, shape=Host), size=3,stroke=0,show.legend = T) +
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
        legend.text = element_text(size=8), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,8))+
  geom_cladelabel(node = 776, label = "bold(Picornaviridae)",parse=T,hjust='center',offset=0.15, fontsize = 2, color="turquoise4") +
  geom_cladelabel(node = 808, label = "bold(Caliciviridae)",parse=T,hjust='center', offset=0.12,fontsize = 2, color="royalblue4") +
  geom_cladelabel(node = 440, label = "bold(Picornaviridae)",parse=T,hjust='center', offset=0.15, fontsize = 2, color="turquoise4") +
  geom_cladelabel(node = 446, label = "bold(Picornaviridae)",parse=T,hjust='center', offset=0.15,fontsize = 2, color="turquoise4") +
  geom_cladelabel(node = 499, label = "bold(Caliciviridae)",parse=T,hjust='center', offset=0.1,fontsize = 2, color="royalblue4") +
  geom_cladelabel(node = 496, label = "bold(Picornaviridae)",parse=T,hjust='center', offset=0.15, fontsize = 2, color="turquoise4") +
  geom_cladelabel(node = 492, label = "bold(Picornaviridae)",parse=T,hjust='center', offset=0.15, fontsize = 2, color="turquoise4") +
  geom_cladelabel(node = 466, label = "bold(Caliciviridae)",parse=T,hjust='center', offset=0.1,fontsize = 2, color="royalblue4") +
  geom_cladelabel(node = 457, label = "bold(Picornaviridae)", parse=T,fontsize = 2,hjust='center', offset=0.15, color="turquoise4") +
  geom_cladelabel(node = 475, label = "bold(Picornaviridae)",parse=T, fontsize = 2,hjust='center', offset=0.15, color="turquoise4") +
  geom_cladelabel(node = 479, label = "bold(Picornaviridae)",parse=T, fontsize = 2,hjust='center', offset=0.15, color="turquoise4") +
  geom_cladelabel(node = 743, label = "bold(Caliciviridae)",parse=T, fontsize = 2,hjust='center', offset=0.1, color="royalblue4") +
  geom_cladelabel(node = 729, label = "bold(Secoviridae/Iflaviridae)",parse=T,hjust='center', offset=0.25,fontsize = 2, color="deeppink4") +
  geom_cladelabel(node = 721, label = "bold(Marnaviridae/Picornaviridae/Iflaviridae)",parse=T,hjust='center', offset=0.45, fontsize = 2, color="deeppink4")+
  geom_cladelabel(node = 533, label = "bold(Marnaviridae/Caliciviridae)",parse=T,hjust='center', offset=0.3, fontsize = 2,color="deeppink4") +
  geom_cladelabel(node = 693, label = "bold(Unclassified)",parse=T, fontsize = 2,hjust='center', offset=0.1, color="deeppink4") +
  geom_cladelabel(node = 668, label = "bold(Unclassified)",parse=T, fontsize = 2,hjust='center', offset=0.12, color="deeppink4") +
  geom_cladelabel(node = 701, label = "bold(Unclassified)",parse=T, fontsize = 2,hjust='center', offset=0.12,color="deeppink4") +
  geom_cladelabel(node = 716, label = "bold(Caliciviridae)",parse=T, fontsize = 2,hjust='center', offset=0.12,color="royalblue4") +
  geom_cladelabel(node = 645, label = "bold(Dicistroviridae/Iflaviridae)",parse=T,hjust='center', offset=0.3, fontsize = 2, color="deeppink4")+
  geom_cladelabel(node = 629, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset=0.15, color="goldenrod4") +
  geom_cladelabel(node = 624, label = "bold(Picornaviridae)",parse=T,hjust='center', offset=0.16, fontsize = 2, color="turquoise4") +
  geom_cladelabel(node = 602, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset=0.15, color="goldenrod4") +
  geom_cladelabel(node = 583, label = "bold(Secoviridae)", parse=T,fontsize = 2,hjust='center', offset=0.15, color="goldenrod4")+
  geom_cladelabel(node = 574, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset=0.15, color="goldenrod4") +
  geom_cladelabel(node = 560, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset=0.15, color="goldenrod4") +
  geom_cladelabel(node = 572, label = "bold(Secoviridae)", parse=T,fontsize = 2,hjust='center', offset=0.15, color="goldenrod4") +
  geom_cladelabel(node = 565, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset=0.15, color="goldenrod4") +
  geom_cladelabel(node = 567, label = "bold(Secoviridae)", parse=T,fontsize = 2,hjust='center', offset=0.17, color="goldenrod4")
p2


#collapse the labeled clades
p3<-collapse(p2, 776)+geom_point2(aes(subset=(node==776)), size=3, shape=22, fill="turquoise1")
p4<-collapse(p3, 808)+geom_point2(aes(subset=(node==808)), size=3, shape=22, fill="royalblue3")
p5<-collapse(p4, 440)+geom_point2(aes(subset=(node==440)), size=3, shape=22, fill="turquoise1")
p6<-collapse(p5, 446)+geom_point2(aes(subset=(node==446)), size=3, shape=22, fill="turquoise1")
p7<-collapse(p6, 499)+geom_point2(aes(subset=(node==499)), size=3, shape=22, fill="royalblue3")
p8<-collapse(p7, 496)+geom_point2(aes(subset=(node==496)), size=3, shape=22, fill="turquoise1")
p9<-collapse(p8, 492)+geom_point2(aes(subset=(node==492)), size=3, shape=22, fill="turquoise1")
p10<-collapse(p9, 466)+geom_point2(aes(subset=(node==466)), size=3, shape=22, fill="royalblue3")
p11<-collapse(p10, 457)+geom_point2(aes(subset=(node==457)), size=3, shape=22, fill="turquoise1")
p12<-collapse(p11, 475)+geom_point2(aes(subset=(node==475)), size=3, shape=22, fill="turquoise1")
p13<-collapse(p12, 479)+geom_point2(aes(subset=(node==479)), size=3, shape=22, fill="turquoise1")
p14<-collapse(p13, 743)+geom_point2(aes(subset=(node==743)), size=3, shape=22, fill="royalblue3")
p15<-collapse(p14, 729)+geom_point2(aes(subset=(node==729)), size=3, shape=22, fill="deeppink4")
p16<-collapse(p15, 721)+geom_point2(aes(subset=(node==721)), size=3, shape=22, fill="deeppink4")
p17<-collapse(p16, 533)+geom_point2(aes(subset=(node==533)), size=3, shape=22, fill="deeppink4")
p18<-collapse(p17, 693)+geom_point2(aes(subset=(node==693)), size=3, shape=22, fill="deeppink4")
p19<-collapse(p18, 668)+geom_point2(aes(subset=(node==668)), size=3, shape=22, fill="deeppink4")
p20<-collapse(p19, 701)+geom_point2(aes(subset=(node==701)), size=3, shape=22, fill="deeppink4")
p21<-collapse(p20, 716)+geom_point2(aes(subset=(node==716)), size=3, shape=22, fill="royalblue3")
p22<-collapse(p21, 645)+geom_point2(aes(subset=(node==645)), size=3, shape=22, fill="deeppink4")
p23<-collapse(p22, 629)+geom_point2(aes(subset=(node==629)), size=3, shape=22, fill="goldenrod1")
p24<-collapse(p23, 624)+geom_point2(aes(subset=(node==624)), size=3, shape=22, fill="turquoise1")
p25<-collapse(p24, 602)+geom_point2(aes(subset=(node==602)), size=3, shape=22, fill="goldenrod1")
p26<-collapse(p25, 583)+geom_point2(aes(subset=(node==583)), size=3, shape=22, fill="goldenrod1")
p27<-collapse(p26, 574)+geom_point2(aes(subset=(node==574)), size=3, shape=22, fill="goldenrod1")
p28<-collapse(p27, 560)+geom_point2(aes(subset=(node==560)), size=3, shape=22, fill="goldenrod1")
p29<-collapse(p28, 572)+geom_point2(aes(subset=(node==572)), size=3, shape=22, fill="goldenrod1")
p30<-collapse(p29, 565)+geom_point2(aes(subset=(node==565)), size=3, shape=22, fill="goldenrod1")
p31<-collapse(p30, 567)+geom_point2(aes(subset=(node==567)), size=3, shape=22, fill="goldenrod1")
p31


##add bootstrap values to this tree
p31.dat <- p31$data
p31.dat$Bootstrap <- NA
Bootstrap<-p31.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p31.dat$label)] <- as.numeric(p31.dat$label[(length(tree.dat$tip_label)+1):length(p0.dat$label)])#fill with label

p31.1 <- p31  %<+% p0.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"))
p31.1






##Clades to collapse
#776 Picornaviridae: Enterovirus with some random iflavirus
#808 Caliciviridae: Norovirus, vesivirus, calicivirus/sadwavirus, dicipivirus, iflavirus
#440 Picornaviridae: cosavirus
#446 picornaviridae: avisivirus, aalivirus, parechovirus, avihepatoivirus mainly avian picornaviridae
#499 caliciviridae: norovirus, calicivirus, nebovirus, one bat sapovirus
#496 unclassified picornaviridae
#492 unclassified picornaviridae
#466 caliciviridae: vesivirus, 3 random viruses
#457 picornaviridae: cardiovirus, malagasivirus, orivirus, tremovirus, random sadwavirus
#475 picornaviridae: limnipivirus, apthovirus, potamipivirus
#479 picornaviridae: megrivirus, gallivirus, oscivirus, rosavirus
#743 caliciviridae: lagosvirus, sapovirus, two random apthovirus
#729: random secoviridae, iflaviridae
#721 random marnaviridae, picornaviridae, iflaviridae
#533 random marnaviridae and caliciviridae
#693 random
#668 random
#701 unclassified
#716 caliciviridae: aparavirus
#645 dicistroviridae and iflaviridae mostly: iflavirus, dicistrovirus, cripavirus, centovirus, triatovirus
#629 secoviridae: nepovirus, and some picornaviridae: megrivirus
#624 unclassified picornaviridae with some iflavirus
#602 mostly secoviridae nepovirus
#583 secoviridae torradovirus
#574 secoviridae: comovirus
#560 secoviridae: nepovirus, one sogarnavirus
#572 secoviridae: waikavirus
#565: secoviridae and polycipiviridae
#567: secoviridae: comovirus




