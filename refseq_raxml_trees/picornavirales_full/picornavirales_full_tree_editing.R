rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)

homewd= "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/"

##Picornavirales all full sequences
setwd(paste0(homewd,"/refseq_raxml_trees/picornavirales_full"))

#load the tree
tree <-  read.tree("T10.raxml.supportFBP") 

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("picornavirales_full_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #408
length(tree$tip.label) #408

#check subgroup names
unique(dat$Family)

colz = c("Caliciviridae" = "royalblue3",    "Picornaviridae"  = "turquoise1",   "Secoviridae"  = "goldenrod1",   "Dicistroviridae"   = "dodgerblue1" ,   "Iflaviridae" = "firebrick1" ,
         "Polycipiviridae"  = "lightpink1" ,    "Marnaviridae"  = "hotpink1" ,  "Solinviviridae" = "lightskyblue" ,   "Coronaviridae"  = "black", "Unclassified"  = "darkorange1")

#pick order for the labels
dat$Family <- factor(dat$Family, levels = c("Dicistroviridae",    "Secoviridae",   "Marnaviridae",   "Picornaviridae",   "Solinviviridae",
                                            "Caliciviridae",    "Iflaviridae",  "Coronaviridae",   "Polycipiviridae",  "Unclassified"))   

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

dat$Host
dat$Host[dat$Host==""] <- NA
dat$Country #some are blank, so convert those to NA
dat$Country[dat$Country==""] <- NA
dat$Collection_Date #these are messy, some are years and some are full dates. I just want years, will manually fix


#now, select the name based on what components are present for each sample

#all components with values:
dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                           dat$Genus[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                           dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|", 
                                                                                                           dat$Host[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                           dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                           dat$Collection_Date[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)])

#and if there is an NA just drop it

#here NA in Isolate only:
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|", 
                                                                                                          dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|", 
                                                                                                          #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|", 
                                                                                                          dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                          dat$Country[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                          dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)])

#and Host only:
dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) &is.na(dat$Host) &!is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)], "|", 
                                                                                                          dat$Genus[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)], "|", 
                                                                                                          dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)], "|", 
                                                                                                          #dat$Host[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                          dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                          dat$Collection_Date[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)])


#and Country only
dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)] <- paste(dat$Accession[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|", 
                                                                                                           dat$Genus[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|", 
                                                                                                           dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|", 
                                                                                                           dat$Host[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|",
                                                                                                           #dat$Country[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                           dat$Collection_Date[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)])

#and accession number only
#check:
dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) &!is.na(dat$Host) & !is.na(dat$Country)] #none with this condition



#now NA in Isolate + Host
#first check if there are any
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)] #no
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)]<- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)], "|",
#                                                                                                        dat$Family[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)], "|",
#                                                                                                        #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|",
#                                                                                                        #dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|",
#                                                                                                        dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)], "|",
#                                                                                                        dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)])
#Isolate + accession only
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)] #no

#Isolate + Country only
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] #no
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
#                                                                                                                 dat$Family[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
#                                                                                                           #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|", 
#                                                                                                           dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|",
#                                                                                                           #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)], "|",
#                                                                                                           dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)])

#now look at the other 2-way NAs: accession + Host, accession + Country, Host+Country
dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$Host) & !is.na(dat$Country)] #none

dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) & is.na(dat$Country)] #none
dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)] #none

#replace
# dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)]<- paste(dat$Accession[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)], "|", 
#                                                                                                                 dat$Family[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)], "|", 
#                                                                                                                 dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)], "|", 
#                                                                                                                 #dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|",
#                                                                                                                 #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)], "|",
#                                                                                                                 dat$Collection_Date[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)])

#and 3-way NAs: Isolate + accession + Host
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$Host) & !is.na(dat$Country)] #none
#Isolate + accession + Country
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) & is.na(dat$Country)] #none
#Isolate Host Country
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)] #none
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)], "|",
#                                                                                                                 dat$Family[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)], "|",
#                                                                                                                 #dat$Isolate[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)], "|",
#                                                                                                                 #dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|",
#                                                                                                                 #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)], "|",
#                                                                                                                 dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)])
# #accession Country Host
dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)] #none
#and all four
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)] #none

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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, bat_Host, Country, Collection_Date, Family, Genus, novel, old_tip_label)

rooted.tree$tip.label <- tree.dat$tip_label

head(tree.dat)
#verify that the old labels match the new

cbind(tree.dat$old_tip_label, rooted.tree$tip.label)

#check out the labels
tree.dat$tip_label#all look good

tree.dat$bat_Host[tree.dat$bat_Host==0] <- "Non-bat host"
tree.dat$bat_Host[tree.dat$bat_Host==1] <- "Bat host"
tree.dat$bat_Host <- as.factor(tree.dat$bat_Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Family, shape=bat_Host), size=3, show.legend = T) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = -.1) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", family="Helvetica", label.size = 0, alpha=.3, size=2, show.legend=FALSE) +#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "left", 
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,8))
p1

##Get the clade numbers so we can collapse unnnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)


#collapsed tree

#add clade labels
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Family, shape=bat_Host), size=3, show.legend = T) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = -.1) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", family="Helvetica", label.size = 0, alpha=.3, size=2, show.legend=FALSE) +#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "left", 
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.2, "cm")) +
  xlim(c(0,8))+
  geom_cladelabel(node = 529, label = "bold(Picornaviridae)",parse=T,hjust='center',offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 565, label = "bold(Caliciviridae)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 569, label = "bold(Picornaviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 580, label = "bold(Caliciviridae)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 618, label = "bold(Picornaviridae)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 610, label = "bold(Picornaviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 524, label = "bold(Picornaviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 494, label = "bold(Picornaviridae/Caliciviridae)",parse=T,hjust='center', offset.text=0.4,fontsize = 2, color="black") +
  geom_cladelabel(node = 652, label = "bold(Calicivirdae)", parse=T,fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 688, label = "bold(Iflaviridae/Polycipiviridae/Picornaviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.6, color="black") +
  geom_cladelabel(node = 460, label = "bold(Iflaviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 483, label = "bold(Dicistroviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 442, label = "bold(Picornaviridae)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 439, label = "bold(Iflaviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black")+
  geom_cladelabel(node = 428, label = "bold(Picornaviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2,color="black") +
  geom_cladelabel(node = 453, label = "bold(Unclassified)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 806, label = "bold(Dicistroviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 804, label = "bold(Iflaviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2,color="black") +
  geom_cladelabel(node = 796, label = "bold(Unclassified)",parse=T, fontsize = 2,hjust='center', offset.text=0.2,color="black") +
  geom_cladelabel(node = 810, label = "bold(Secoviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black")+
  geom_cladelabel(node = 698, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 789, label = "bold(Iflaviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 734, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 741, label = "bold(Unclassified)", parse=T,fontsize = 2,hjust='center', offset.text=0.2, color="black")+
  geom_cladelabel(node = 743, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 762, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 784, label = "bold(Secoviridae)", parse=T,fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 776, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black")
p2


#collapse the labeled clades
p3<-collapse(p2, 529)+geom_point2(aes(subset=(node==529)), size=3, shape=15, color="turquoise1", alpha=0.99)
p4<-collapse(p3, 565)+geom_point2(aes(subset=(node==565)), size=3, shape=15, color="royalblue3")
p5<-collapse(p4, 569)+geom_point2(aes(subset=(node==569)), size=3, shape=15, color="turquoise1")
p6<-collapse(p5, 580)+geom_point2(aes(subset=(node==580)), size=3, shape=15, color="royalblue3")
p7<-collapse(p6, 618)+geom_point2(aes(subset=(node==618)), size=3, shape=15, color="turquoise1")
p8<-collapse(p7, 610)+geom_point2(aes(subset=(node==610)), size=3, shape=15, color="turquoise1")
p9<-collapse(p8, 524)+geom_point2(aes(subset=(node==524)), size=3, shape=15, color="turquoise1")
p10<-collapse(p9, 494)+geom_point2(aes(subset=(node==494)), size=3, shape=15, color="deeppink4")
p11<-collapse(p10, 652)+geom_point2(aes(subset=(node==652)), size=3, shape=15, color="royalblue3")
p12<-collapse(p11, 688)+geom_point2(aes(subset=(node==688)), size=3, shape=15, color="deeppink4")
p13<-collapse(p12, 460)+geom_point2(aes(subset=(node==460)), size=3, shape=15, color="firebrick1")
p14<-collapse(p13, 483)+geom_point2(aes(subset=(node==483)), size=3, shape=15, color="dodgerblue1")
p15<-collapse(p14, 442)+geom_point2(aes(subset=(node==442)), size=3, shape=15, color="deeppink4")
p16<-collapse(p15, 439)+geom_point2(aes(subset=(node==439)), size=3, shape=15, color="firebrick1")
p17<-collapse(p16, 428)+geom_point2(aes(subset=(node==428)), size=3, shape=15, color="deeppink4")
p18<-collapse(p17, 453)+geom_point2(aes(subset=(node==453)), size=3, shape=15, color="deeppink4")
p19<-collapse(p18, 806)+geom_point2(aes(subset=(node==806)), size=3, shape=15, color="dodgerblue1")
p20<-collapse(p19, 804)+geom_point2(aes(subset=(node==804)), size=3, shape=15, color="firebrick1")
p21<-collapse(p20, 796)+geom_point2(aes(subset=(node==796)), size=3, shape=15, color="deeppink4")
p22<-collapse(p21, 810)+geom_point2(aes(subset=(node==810)), size=3, shape=15, color="deeppink4")
p23<-collapse(p22, 698)+geom_point2(aes(subset=(node==698)), size=3, shape=15, color="goldenrod1")
p24<-collapse(p23, 789)+geom_point2(aes(subset=(node==789)), size=3, shape=15, color="firebrick1")
p25<-collapse(p24, 734)+geom_point2(aes(subset=(node==734)), size=3, shape=15, color="goldenrod1")
p26<-collapse(p25, 741)+geom_point2(aes(subset=(node==741)), size=3, shape=15, color="deeppink4")
p27<-collapse(p26, 743)+geom_point2(aes(subset=(node==743)), size=3, shape=15, color="goldenrod1")
p28<-collapse(p27, 762)+geom_point2(aes(subset=(node==762)), size=3, shape=15, color="goldenrod1")
p29<-collapse(p28, 784)+geom_point2(aes(subset=(node==784)), size=3, shape=15, color="goldenrod1")
p30<-collapse(p29, 776)+geom_point2(aes(subset=(node==776)), size=3, shape=15, color="goldenrod1")
p30


##Clades to collapse
#529 Picornaviridae
#565 caliciviridae
#569 picornaviridae
#580 caliciviridae
#618 picornaviridae
#610 picornaviridae
#524 picornaviridae
#494 picornaviridae/caliciviridae/random iflaviridae
#652 caliciviridae
#688 iflaviridae/polycipiviridae/picornaviridae
#460 iflaviridae
#483 dicistroviridae
#442 picornaviridae
#439 iflaviridae
#428 picornaviridae
#453 unclassified
#806 dicistroviridae
#804 iflaviridae
#796 Unclassified
#810 secoviridae
#698 secoviridae
#789 iflaviridae
#734 secoviridae
#741 unclassified
#743 secoviridae
#762 secoviridae
#784 secoviridae
#776 secoviridae








