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
setwd(paste0(homewd,"/raxml_trees/secoviridae"))

#load the tree
tree <-  read.tree("T1.raxml.supportFBP") 

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("secoviridae_refseq_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #84
length(tree$tip.label) #84

#check subgroup names
unique(dat$Genus)

colz = c("Waikavirus" = "royalblue3",    "Fabavirus"  = "turquoise1",   "Sadwavirus"  = "goldenrod1",   "Comovirus"   = "dodgerblue1" ,   "Nepovirus" = "firebrick1" ,
         "Sequivirus"  = "lightpink1" ,    "Cheravirus"  = "hotpink1" ,  "Stralarivirus" = "lightskyblue" ,   "Coronavirus"  = "black", "Torradovirus"  = "darkorange1")

#pick order for the labels
dat$Genus <- factor(dat$Genus, levels = c("Waikavirus" ,  "Fabavirus",   "Sadwavirus",   "Comovirus",   "Nepovirus",
                                          "Sequivirus",    "Cheravirus",  "Stralarivirus",  "Torradovirus",  "Coronavirus"))   

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
#                                                                                                        dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)], "|",
#                                                                                                        #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|",
#                                                                                                        #dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|",
#                                                                                                        dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)], "|",
#                                                                                                        dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) &!is.na(dat$Country)])
#Isolate + accession only
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)] #no

#Isolate + Country only
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] #no
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
#                                                                                                                 dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
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
#                                                                                                                 dat$Genus[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)], "|", 
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
#                                                                                                                 dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$Host) & is.na(dat$Country)], "|",
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
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, bat_Host, Country, Collection_Date, Genus, Genus, novel, old_tip_label)

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
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=bat_Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Genus="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=2,  nudge_x=0.1) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "left", 
        legend.direction = "vertical",
        legend.text = element_text(size=8), 
        legend.key.size = unit(0.3, "cm")) +
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
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=bat_Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Genus="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=2,  nudge_x=0.1) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=8), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,8))+
  geom_cladelabel(node = 132, label = "bold(Nepovirus)",parse=T,hjust='center',offset=0.3, fontsize = 2, color="firebrick3") +
  geom_cladelabel(node = 114, label = "bold(Sadwavirus)",parse=T,hjust='center', offset=0.3,fontsize = 2, color="goldenrod4") +
  geom_cladelabel(node = 160, label = "bold(Torradovirus)",parse=T,hjust='center', offset=0.3, fontsize = 2, color="darkorange3") +
  geom_cladelabel(node = 156, label = "bold(Cheravirus)",parse=T,hjust='center', offset=0.3,fontsize = 2, color="hotpink3") +
  geom_cladelabel(node = 154, label = "bold(Sequivirus)",parse=T,hjust='center', offset=0.3,fontsize = 2, color="lightpink3") +
  geom_cladelabel(node = 87, label = "bold(Fabavirus)",parse=T,hjust='center', offset=0.3, fontsize = 2, color="turquoise4") +
  geom_cladelabel(node = 94, label = "bold(Comovirus)",parse=T,hjust='center', offset=0.3, fontsize = 2, color="dodgerblue3")
p2


#collapse the labeled clades
p3<-collapse(p2, 132)+geom_point2(aes(subset=(node==132)), size=3, shape=22, fill="firebrick1")
p4<-collapse(p3, 114)+geom_point2(aes(subset=(node==114)), size=3, shape=22, fill="goldenrod1")
p5<-collapse(p4, 160)+geom_point2(aes(subset=(node==160)), size=3, shape=22, fill="darkorange1")
p6<-collapse(p5, 156)+geom_point2(aes(subset=(node==156)), size=3, shape=22, fill="hotpink1")
p7<-collapse(p6, 154)+geom_point2(aes(subset=(node==154)), size=3, shape=22, fill="lightpink1")
p8<-collapse(p7, 87)+geom_point2(aes(subset=(node==87)), size=3, shape=22, fill="turquoise1")
p9<-collapse(p8, 94)+geom_point2(aes(subset=(node==94)), size=3, shape=22, fill="dodgerblue1")
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

##Clades to collapse
##132 Nepovirus
#114 sadwavirus
#160 Torradovirus
#156 cheravirus
#154 sequivirus
#87 fabavirus
#94 Comovirus