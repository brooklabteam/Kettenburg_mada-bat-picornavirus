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
setwd(paste0(homewd,"/raxml_trees/picornaviridae"))

#load the tree
tree <-  read.tree("T10.raxml.supportFBP") 

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("picornaviridae_refseq_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #179
length(tree$tip.label) #179

#check subgroup names
unique(dat$Genus)

# colz = c("Cardiovirus" = "cadetblue",    "Enterovirus"  = "cadetblue1",   "Hepatovirus"  = "cadetblue2",   
#          "Kobuvirus"   = "cadetblue3" ,   "Parechovirus" = "chartreuse" ,"Erbovirus"  = "chartreuse3" ,    
#          "Teschovirus"  = "chartreuse4" ,  "Sapelovirus" = "chocolate" ,   "Tremovirus"  = "chocolate2" ,   
#          "Anativirus"   = "chocolate3" ,"Avihepatovirus"  = "chocolate4","Aquamavirus"  = "coral",   
#          "Aphthovirus"  = "coral2" ,  "Senecavirus"  = "coral3" ,  "Cosavirus"  = "cyan",
#          "Salivirus"   = "cyan3",    "Passerivirus"  = "cyan4", "Oscivirus"   = "darkgoldenrod1" , 
#          "Unclassified picornavirus"= "thistle1"  ,   "Mischivirus" = "darkgoldenrod3","Pasivirus"   = "firebrick1"  ,  
#          "Gallivirus"   = "firebrick3" ,  "Limnipivirus"  = "firebrick4",  "Hunnivirus"   = "darkolivegreen1",   
#          "Dicipivirus" = "darkolivegreen3","Megrivirus"  = "darkorchid1" ,   "Potamipivirus"  = "darkorchid4", 
#          "Sakobuvirus"  = "darkseagreen1" ,  "Sicinivirus"  = "darkseagreen3" ,  "Aalivirus"  = "deeppink",
#          "Mosavirus"  = "deeppink3" ,    "Rosavirus"  = "deeppink4" ,    "Avisivirus"   = "hotpink",   
#          "Orivirus"   = "hotpink3" ,    "Crohivirus" = "hotpink4","Torchivirus" = "indianred1" ,   
#          "Bopivirus"  = "indianred3"  ,   "Malagasivirus" = "indianred4" , "Harkavirus"  = "mediumpurple1" ,   
#          "Ampivirus"  = "mediumpurple3" ,"Livupivirus" = "mediumpurple4" ,   "Kunsagivirus"  = "orange",  
#          "Shanbavirus"  = "orange3" ,  "Rafivirus"   = "seagreen1",  "Coronavirus" ="black",  
#          "Poecivirus"  = "seagreen3" ,"Rabovirus"   = "seagreen4",    "Tottorivirus"  = "skyblue1" , 
#          "Ailurivirus" = "skyblue3", "Madagascar bat kobuvirus" ="cadetblue3", "Bat picornavirus"="yellow", "Picornavirus"="grey2")


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

#and accession number only
#check:
# dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) &!is.na(dat$source) & !is.na(dat$Country)] #none with this condition
# 
# 
# 
# #now NA in Isolate + source
# #first check if there are any
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)] #no
# # dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)]<- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
# #                                                                                                        dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
# #                                                                                                        #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|",
# #                                                                                                        #dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|",
# #                                                                                                        dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)], "|",
# #                                                                                                        dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) &!is.na(dat$Country)])
# #Isolate + accession only
# dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) &!is.na(dat$Country)] #no
# 
# #Isolate + Country only
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)] #no
# # dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|", 
# #                                                                                                                 dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|", 
# #                                                                                                           #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$source) & is.na(dat$Country)], "|", 
# #                                                                                                           dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|",
# #                                                                                                           #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) &!is.na(dat$Country)], "|",
# #                                                                                                           dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)])
# 
# #now look at the other 2-way NAs: accession + source, accession + Country, source+Country
# dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$source) & !is.na(dat$Country)] #none
# 
# dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) & is.na(dat$Country)] #none
# dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)] #none
# 
# #replace
# # dat$new_label[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)]<- paste(dat$Accession[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|", 
# #                                                                                                                 dat$Genus[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|", 
# #                                                                                                                 dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|", 
# #                                                                                                                 #dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|",
# #                                                                                                                 #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) &!is.na(dat$Country)], "|",
# #                                                                                                                 dat$Collection_Date[!is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)])
# 
# #and 3-way NAs: Isolate + accession + source
# dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$source) & !is.na(dat$Country)] #none
# #Isolate + accession + Country
# dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) & is.na(dat$Country)] #none
# #Isolate source Country
# dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)] #none
# # dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|",
# #                                                                                                                 dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|",
# #                                                                                                                 #dat$Isolate[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)], "|",
# #                                                                                                                 #dat$source[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$source) &is.na(dat$Country)], "|",
# #                                                                                                                 #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$source) &!is.na(dat$Country)], "|",
# #                                                                                                                 dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)])
# # #accession Country source
# dat$new_label[!is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)] #none
# #and all four
# dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & is.na(dat$source) & is.na(dat$Country)] #none

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
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Genus="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=2, nudge_x=0.05) +
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
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Genus="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=2, nudge_x=0.05) +
  guides(fill="none")+#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=0,y=-3, linesize = .5) +
  theme(legend.position = "none", 
        legend.direction = "vertical",
        legend.text = element_text(size=8), 
        legend.key.size = unit(0.3, "cm")) +
  xlim(c(0,8))+
  geom_cladelabel(node = 261, label = "bold(Enterovirus)",parse=T,hjust='center',offset=0.3, fontsize = 2, color="cadetblue4") +
  geom_cladelabel(node = 235, label = "bold(Hepatovirus)",parse=T,hjust='center', offset=0.3,fontsize = 2, color="cadetblue4") +
  geom_cladelabel(node = 213, label = "bold(Avisivirus)",parse=T,hjust='center', offset=0.3, fontsize = 2, color="hotpink2") +
  geom_cladelabel(node = 216, label = "bold(Limnipivirus/Potamipivirus)",parse=T,hjust='center', offset=0.55,fontsize = 2, color="deeppink4") +
  geom_cladelabel(node = 205, label = "bold(Parechovirus/Shanbavirus/Avihepatovirus)",parse=T,hjust='center', offset=0.8,fontsize = 2, color="deeppink4") +
  geom_cladelabel(node = 336, label = "bold(Ampivirus)",parse=T,hjust='center', offset=0.3, fontsize = 2, color="pink2") +
  geom_cladelabel(node = 190, label = "bold(Malagasivirus/Tottorivirus/Hunnivirus)",parse=T,hjust='center', offset=0.7, fontsize = 2, color="deeppink4") +
  geom_cladelabel(node = 188, label = "bold(Mosavirus)",parse=T,hjust='center', offset=0.3,fontsize = 2, color="lightblue4") +
  geom_cladelabel(node = 185, label = "bold(Cosavirus)", parse=T,fontsize = 2,hjust='center', offset=0.3, color="darkgoldenrod4") +
  geom_cladelabel(node = 342, label = "bold(Megrivirus/Ailurivirus)",parse=T, fontsize = 2,hjust='center', offset=0.45, color="deeppink4") +
  geom_cladelabel(node = 346, label = "bold(Cardiovirus)",parse=T, fontsize = 2,hjust='center', offset=0.3, color="4") +
  geom_cladelabel(node = 354, label = "bold(Cardiovirus/Senecavirus)",parse=T, fontsize = 2,hjust='center', offset=0.5, color="deeppink4") +
  geom_cladelabel(node = 351, label = "bold(Apthovirus)",parse=T,hjust='center', offset=0.27,fontsize = 2, color="darkgoldenrod2") +
  geom_cladelabel(node = 353, label = "bold(Bopivirus/Erbovirus)",parse=T,hjust='center', offset=0.4, fontsize = 2, color="deeppink4")+
  geom_cladelabel(node = 333, label = "bold(Rosavirus)",parse=T,hjust='center', offset=0.27, fontsize = 2,color="hotpink1") +
  geom_cladelabel(node = 332, label = "bold(Harkavirus/Dicipivirus)",parse=T, fontsize = 2,hjust='center', offset=0.45, color="deeppink4") +
  geom_cladelabel(node = 326, label = "bold(Oscivirus/Livupivirus/Rafivirus)",parse=T, fontsize = 2,hjust='center', offset=0.6, color="deeppink4") +
  geom_cladelabel(node = 323, label = "bold(Megrivirus/Parechovirus)",parse=T, fontsize = 2,hjust='center', offset=0.5,color="deeppink4") +
  geom_cladelabel(node = 295, label = "bold(Gallivirus/Megrivirus)",parse=T, fontsize = 2,hjust='center', offset=0.4,color="deeppink4") +
  geom_cladelabel(node = 303, label = "bold(Passerivirus/Sicinvirus)",parse=T,hjust='center', offset=0.5, fontsize = 2, color="deeppink4")+
  geom_cladelabel(node = 320, label = "bold(Salivirus/Sakobuvirus)",parse=T, fontsize = 2,hjust='center', offset=0.4, color="deeppink4")
p2


#collapse the labeled clades
p3<-collapse(p2, 261)+geom_point2(aes(subset=(node==261)), size=3, shape=22, fill="cadetblue2")
p4<-collapse(p3, 235)+geom_point2(aes(subset=(node==235)), size=3, shape=22, fill="cadetblue3")
p5<-collapse(p4, 213)+geom_point2(aes(subset=(node==213)), size=3, shape=22, fill="hotpink2")
p6<-collapse(p5, 216)+geom_point2(aes(subset=(node==216)), size=3, shape=22, fill="deeppink4")
p7<-collapse(p6, 205)+geom_point2(aes(subset=(node==205)), size=3, shape=22, fill="deeppink4")
p8<-collapse(p7, 336)+geom_point2(aes(subset=(node==336)), size=3, shape=22, fill="pink2")
p9<-collapse(p8, 190)+geom_point2(aes(subset=(node==190)), size=3, shape=22, fill="deeppink4")
p10<-collapse(p9, 188)+geom_point2(aes(subset=(node==188)), size=3, shape=22, fill="lightblue4")
p11<-collapse(p10, 185)+geom_point2(aes(subset=(node==185)), size=3, shape=22, fill="darkgoldenrod4")
p12<-collapse(p11, 342)+geom_point2(aes(subset=(node==342)), size=3, shape=22, fill="deeppink4")
p13<-collapse(p12, 346)+geom_point2(aes(subset=(node==346)), size=3, shape=22, fill="cadetblue1")
p14<-collapse(p13, 354)+geom_point2(aes(subset=(node==354)), size=3, shape=22, fill="deeppink4")
p15<-collapse(p14, 351)+geom_point2(aes(subset=(node==351)), size=3, shape=22, fill="darkgoldenrod2")
p16<-collapse(p15, 353)+geom_point2(aes(subset=(node==353)), size=3, shape=22, fill="deeppink4")
p17<-collapse(p16, 333)+geom_point2(aes(subset=(node==333)), size=3, shape=22, fill="hotpink1")
p18<-collapse(p17, 332)+geom_point2(aes(subset=(node==332)), size=3, shape=22, fill="deeppink4")
p19<-collapse(p18, 326)+geom_point2(aes(subset=(node==326)), size=3, shape=22, fill="deeppink4")
p20<-collapse(p19, 323)+geom_point2(aes(subset=(node==323)), size=3, shape=22, fill="deeppink4")
p21<-collapse(p20, 295)+geom_point2(aes(subset=(node==295)), size=3, shape=22, fill="deeppink4")
p22<-collapse(p21, 303)+geom_point2(aes(subset=(node==303)), size=3, shape=22, fill="deeppink4")
p23<-collapse(p22, 320)+geom_point2(aes(subset=(node==320)), size=3, shape=22, fill="deeppink4")
p23


##add bootstrap values to this tree
p23.dat <- p23$data
p23.dat$Bootstrap <- NA
Bootstrap<-p23.dat$Bootstrap[(length(tree.dat$tip_label)+1):length(p23.dat$label)] <- as.numeric(p23.dat$label[(length(tree.dat$tip_label)+1):length(p23.dat$label)])#fill with label

p23.1 <- p23  %<+% p23.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=Bootstrap, show.legend = T), shape=21, stroke=0)+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100))+
  #guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1))+
  theme(legend.position = "left",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"))
p23.1






##Clades to collapse
#261 enterovirus
#235 hepatovirus
#213 avisivirus
#216 limnipivirus/potamipivirus
#205 parechovirus/shanbavirus/avihepatovirus
#336 Ampivirus
#190 malagasivirus/tottorivirus/hunnivirus
#188 mosavirus
#185 cosavirus
#342 megrivirus/ailurivirus
#346 cardiovirus
#354 cardiovirus/senecavirus
#351 apthovirus
#353 bopivirus/erbovirus
#333 rosavirus
#332 harkavirus/dicipivirus
#326 oscivirus/livupivirus/rafivirus
#323 megrivirus/parechovirus
#295 gallivirus/megrivirus
#303 passerivirus/sicinvirus
#320 salivirus/sakobuvirus





