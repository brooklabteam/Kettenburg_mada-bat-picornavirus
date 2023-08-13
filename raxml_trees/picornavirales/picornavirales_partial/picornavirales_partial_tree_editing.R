rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)

homewd= "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/"

##Picornavirales all partial sequences
setwd(paste0(homewd,"/refseq_raxml_trees/picornavirales_partial"))

#load the tree
tree <-  read.tree("T10.raxml.supportFBP") 

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("picornavirales_partial_metadata.csv"), header = T, stringsAsFactors = F)
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
  xlim(c(0,7))+
  geom_cladelabel(node = 536, label = "bold(Picornaviridae)",parse=T,hjust='center',offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 641, label = "bold(Secoviridae)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 644, label = "bold(Dicistroviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 497, label = "bold(Picornaviridae)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 486, label = "bold(Dicistroviridae)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 516, label = "bold(Picornaviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 519, label = "bold(Picornaviridae/Dicistroviridae)",parse=T,hjust='center', offset.text=0.4, fontsize = 2, color="black") +
  geom_cladelabel(node = 666, label = "bold(Picornaviridae/Dicistroviridae)",parse=T,hjust='center', offset.text=0.4,fontsize = 2, color="black") +
  geom_cladelabel(node = 669, label = "bold(Picornaviridae/Dicistroviridae)", parse=T,fontsize = 2,hjust='center', offset.text=0.4, color="black") +
  geom_cladelabel(node = 467, label = "bold(Picornaviridae/Dicistroviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.4, color="black") +
  geom_cladelabel(node = 455, label = "bold(Iflaviridae/Picornaviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.4, color="black") +
  geom_cladelabel(node = 771, label = "bold(Dicistroviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 767, label = "bold(Unclassified)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 444, label = "bold(Iflaviridae/Picornaviridae)",parse=T,hjust='center', offset.text=0.4, fontsize = 2, color="black")+
  geom_cladelabel(node = 448, label = "bold(Secoviridae/Iflaviridae)",parse=T,hjust='center', offset.text=0.4, fontsize = 2,color="black") +
  geom_cladelabel(node = 673, label = "bold(Iflaviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 752, label = "bold(Secoviridae/Picornaviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.4, color="black") +
  geom_cladelabel(node = 739, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2,color="black") +
  geom_cladelabel(node = 723, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2,color="black") +
  geom_cladelabel(node = 714, label = "bold(Secoviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black")+
  geom_cladelabel(node = 681, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 710, label = "bold(Secoviridae)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 702, label = "bold(Secoviridae)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 689, label = "bold(Secoviridae)", parse=T,fontsize = 2,hjust='center', offset.text=0.2, color="black")
p2


#collapse the labeled clades
p3<-collapse(p2, 536)+geom_point2(aes(subset=(node==536)), size=3, shape=15, color="turquoise1", alpha=0.99)
p4<-collapse(p3, 641)+geom_point2(aes(subset=(node==641)), size=3, shape=15, color="goldenrod1")
p5<-collapse(p4, 644)+geom_point2(aes(subset=(node==644)), size=3, shape=15, color="royalblue3")
p6<-collapse(p5, 497)+geom_point2(aes(subset=(node==497)), size=3, shape=15, color="turquoise1")
p7<-collapse(p6, 486)+geom_point2(aes(subset=(node==486)), size=3, shape=15, color="royalblue3")
p8<-collapse(p7, 516)+geom_point2(aes(subset=(node==516)), size=3, shape=15, color="turquoise1")
p9<-collapse(p8, 519)+geom_point2(aes(subset=(node==519)), size=3, shape=15, color="deeppink4")
p10<-collapse(p9, 666)+geom_point2(aes(subset=(node==666)), size=3, shape=15, color="deeppink4")
p11<-collapse(p10, 669)+geom_point2(aes(subset=(node==669)), size=3, shape=15, color="deeppink4")
p12<-collapse(p11, 467)+geom_point2(aes(subset=(node==467)), size=3, shape=15, color="deeppink4")
p13<-collapse(p12, 455)+geom_point2(aes(subset=(node==455)), size=3, shape=15, color="deeppink4")
p14<-collapse(p13, 771)+geom_point2(aes(subset=(node==771)), size=3, shape=15, color="dodgerblue1")
p15<-collapse(p14, 767)+geom_point2(aes(subset=(node==767)), size=3, shape=15, color="deeppink4")
p16<-collapse(p15, 444)+geom_point2(aes(subset=(node==444)), size=3, shape=15, color="deeppink4")
p17<-collapse(p16, 448)+geom_point2(aes(subset=(node==448)), size=3, shape=15, color="deeppink4")
p18<-collapse(p17, 673)+geom_point2(aes(subset=(node==673)), size=3, shape=15, color="firebrick1")
p19<-collapse(p18, 752)+geom_point2(aes(subset=(node==752)), size=3, shape=15, color="deeppink4")
p20<-collapse(p19, 739)+geom_point2(aes(subset=(node==739)), size=3, shape=15, color="goldenrod1")
p21<-collapse(p20, 723)+geom_point2(aes(subset=(node==723)), size=3, shape=15, color="goldenrod1")
p22<-collapse(p21, 714)+geom_point2(aes(subset=(node==714)), size=3, shape=15, color="goldenrod1")
p23<-collapse(p22, 681)+geom_point2(aes(subset=(node==681)), size=3, shape=15, color="goldenrod1")
p24<-collapse(p23, 710)+geom_point2(aes(subset=(node==710)), size=3, shape=15, color="goldenrod1")
p25<-collapse(p24, 702)+geom_point2(aes(subset=(node==702)), size=3, shape=15, color="goldenrod1")
p26<-collapse(p25, 689)+geom_point2(aes(subset=(node==689)), size=3, shape=15, color="goldenrod1")
p26


##Clades to collapse
#536 picornaviridae
#641 secoviridae
#644 dicistroviridae
#497 picornaviridae
#486 dicistroviridae
#516 picornaviridae
#519 picornaviridae/dicistroviridae
#666 picornaviridae/dicistroviridae
#669 picornaviridae/dicistroviridae
#467 picornaviridae/dicistroviridae
#455 iflaviridae/picornaviridae
#771 dicistroviridae
#767 Unclassified
#444 iflaviridae/picornaviridae
#448 secoviridae/iflaviridae
#673 iflaviridae
#752 secoviridae/picornaviridae
#739 secoviridae
#723 secoviridae
#714 secoviridae
#681 secoviridae
#710 secoviridae
#702 secoviridae
#689 secoviridae









