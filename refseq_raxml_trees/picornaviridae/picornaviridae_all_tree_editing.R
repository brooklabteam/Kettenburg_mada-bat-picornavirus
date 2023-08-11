rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)

homewd= "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/"


##Picornaviridae all full and partial sequences
setwd(paste0(homewd,"/raxml_trees/picornaviridae_all"))

#load the tree
tree <-  read.tree("T2.raxml.supportFBP")

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("picornaviridae_all_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #181 
length(rooted.tree$tip.label) #181

#check subgroup names
unique(dat$Genus)

colz = c("Cardiovirus" = "cadetblue",    "Enterovirus"  = "cadetblue1",   "Hepatovirus"  = "cadetblue2",   
         "Kobuvirus"   = "cadetblue3" ,   "Parechovirus" = "chartreuse" ,"Erbovirus"  = "chartreuse3" ,    
         "Teschovirus"  = "chartreuse4" ,  "Sapelovirus" = "chocolate" ,   "Tremovirus"  = "chocolate2" ,   
         "Anativirus"   = "chocolate3" ,"Avihepatovirus"  = "chocolate4","Aquamavirus"  = "coral",   
         "Aphthovirus"  = "coral2" ,  "Senecavirus"  = "coral3" ,  "Cosavirus"  = "cyan",
         "Salivirus"   = "cyan3",    "Passerivirus"  = "cyan4", "Oscivirus"   = "darkgoldenrod1" , 
         "Unclassified picornavirus"= "thistle1"  ,   "Mischivirus" = "darkgoldenrod3","Pasivirus"   = "firebrick1"  ,  
         "Gallivirus"   = "firebrick3" ,  "Limnipivirus"  = "firebrick4",  "Hunnivirus"   = "darkolivegreen1",   
         "Dicipivirus" = "darkolivegreen3","Megrivirus"  = "darkorchid1" ,   "Potamipivirus"  = "darkorchid4", 
         "Sakobuvirus"  = "darkseagreen1" ,  "Sicinivirus"  = "darkseagreen3" ,  "Aalivirus"  = "deeppink",
         "Mosavirus"  = "deeppink3" ,    "Rosavirus"  = "deeppink4" ,    "Avisivirus"   = "hotpink",   
         "Orivirus"   = "hotpink3" ,    "Crohivirus" = "hotpink4","Torchivirus" = "indianred1" ,   
         "Bopivirus"  = "indianred3"  ,   "Malagasivirus" = "indianred4" , "Harkavirus"  = "mediumpurple1" ,   
         "Ampivirus"  = "mediumpurple3" ,"Livupivirus" = "mediumpurple4" ,   "Kunsagivirus"  = "orange",  
         "Shanbavirus"  = "orange3" ,  "Rafivirus"   = "seagreen1",  "Coronavirus" ="black",  
         "Poecivirus"  = "seagreen3" ,"Rabovirus"   = "seagreen4",    "Tottorivirus"  = "skyblue1" , 
         "Ailurivirus" = "skyblue3", "Madagascar bat kobuvirus" ="cadetblue3", "Bat picornavirus"="yellow", "Picornavirus"="grey2")

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
                                          "Shanbavirus",  "Rafivirus",  "Coronavirus",  
                                          "Poecivirus", "Rabovirus",    "Tottorivirus", 
                                          "Ailurivirus", "Madagascar bat kobuvirus", "Bat picornavirus", "Picornavirus"))   


dat$novel <- as.factor(dat$novel)

# #rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)
# dat.2 <- cbind.data.frame(tip_label = rooted.tree$tip.label)
# dat.plot <- merge(dat.2, dat, by ="tip_label", all.x = T, sort=F)
#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = "none", legend.title = element_blank())
#quartz(p)
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now check that you don't have NAs and blanks throughout the component parts of the name
#you also can edit the original csv file to replace these blanks if you can find the correct info

dat$Isolate  #some are blank, so convert those to NA
dat$Isolate[dat$Isolate==""] <- NA
dat$Accession #all good

#give false numbers to mada samples
#dat$Accession[is.na(dat$Accession)] <- paste(sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",1), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",2), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",3), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",4), sep="_")

dat$Host #no blanks because I manually did them
#dat$Host[dat$Host==""] <- NA
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

#Isolate + accession only
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)] #no

#Isolate + Country only
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] #yes
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
                                                                                                                dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
                                                                                                          #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|", 
                                                                                                          dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|",
                                                                                                          #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                          dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)])

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
#accession Country Host
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
tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = T, sort = F)

names(tree.dat)

tree.dat$tip_label <- tree.dat$new_label
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, bat_Host, Country, Collection_Date, Genus, novel, old_tip_label)

rooted.tree$tip.label <- tree.dat$tip_label

head(tree.dat)
#verify that the old labels match the new
cbind(tree.dat$old_tip_label, rooted.tree$tip.label) #they match

#check out the labels
tree.dat$tip_label#all look good

tree.dat$bat_Host[tree.dat$bat_Host==0] <- "non-bat Host"
tree.dat$bat_Host[tree.dat$bat_Host==1] <- "bat Host"
tree.dat$bat_Host <- as.factor(tree.dat$bat_Host)
shapez = c("bat Host" =  17, "non-bat Host" = 15)
colz2 = c('1' =  "yellow", '0' = "white")

 
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=bat_Host), size=3, show.legend = T) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = -.08) +
  scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", family="Helvetica", label.size = 0, alpha=.3, size=2, show.legend=F) +#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=3,y=-5, linesize = .5) +
  theme(legend.position = "left", legend.title = element_blank()) +
  xlim(c(0,4))
p1

##collapse unnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)
#200 enterovirus
#350 enterovirus
#242 avihepatovirus, aalivirus, parechovirus, crohivirus, limnipivirus, aquamavirus, pasivirus, avisivirus, orivirus, shanbavirus
#295 rosavirus
#299 apthovirus, erbovirus, bopivirus
#273 cosavirus
#266 cardiovirus
#271 mischivirus
#310 megrivirus
#315 gallivirus
#332 salivirus










##Picornaviridae full and partial p1 region
setwd(paste0(homewd,"/raxml_trees/picornaviridae_all_p1_region"))

#load the tree
tree <-  read.tree("T2.raxml.supportFBP")

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("picornaviridae_p1_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #177
length(rooted.tree$tip.label) #177

#check subgroup names
unique(dat$Genus)

colz = c("Cardiovirus" = "cadetblue",    "Enterovirus"  = "cadetblue1",   "Hepatovirus"  = "cadetblue2",   
         "Kobuvirus"   = "cadetblue3" ,   "Parechovirus" = "chartreuse" ,"Erbovirus"  = "chartreuse3" ,    
         "Teschovirus"  = "chartreuse4" ,  "Sapelovirus" = "chocolate" ,   "Tremovirus"  = "chocolate2" ,   
         "Anativirus"   = "chocolate3" ,"Avihepatovirus"  = "chocolate4","Aquamavirus"  = "coral",   
         "Aphthovirus"  = "coral2" ,  "Senecavirus"  = "coral3" ,  "Cosavirus"  = "cyan",
         "Salivirus"   = "cyan3",    "Passerivirus"  = "cyan4", "Oscivirus"   = "darkgoldenrod1" , 
         "Unclassified picornavirus"= "thistle1"  ,   "Mischivirus" = "darkgoldenrod3","Pasivirus"   = "firebrick1"  ,  
         "Gallivirus"   = "firebrick3" ,  "Limnipivirus"  = "firebrick4",  "Hunnivirus"   = "darkolivegreen1",   
         "Dicipivirus" = "darkolivegreen3","Megrivirus"  = "darkorchid1" ,   "Potamipivirus"  = "darkorchid4", 
         "Sakobuvirus"  = "darkseagreen1" ,  "Sicinivirus"  = "darkseagreen3" ,  "Aalivirus"  = "deeppink",
         "Mosavirus"  = "deeppink3" ,    "Rosavirus"  = "deeppink4" ,    "Avisivirus"   = "hotpink",   
         "Orivirus"   = "hotpink3" ,    "Crohivirus" = "hotpink4","Torchivirus" = "indianred1" ,   
         "Bopivirus"  = "indianred3"  ,   "Malagasivirus" = "indianred4" , "Harkavirus"  = "mediumpurple1" ,   
         "Ampivirus"  = "mediumpurple3" ,"Livupivirus" = "mediumpurple4" ,   "Kunsagivirus"  = "orange",  
         "Shanbavirus"  = "orange3" ,  "Rafivirus"   = "seagreen1",  "Coronavirus" ="black",  
         "Poecivirus"  = "seagreen3" ,"Rabovirus"   = "seagreen4",    "Tottorivirus"  = "skyblue1" , 
         "Ailurivirus" = "skyblue3", "Madagascar bat kobuvirus" ="cadetblue3", "Bat picornavirus"="yellow", "Picornavirus"="grey2")

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
                                          "Shanbavirus",  "Rafivirus",  "Coronavirus",  
                                          "Poecivirus", "Rabovirus",    "Tottorivirus", 
                                          "Ailurivirus", "Madagascar bat kobuvirus", "Bat picornavirus", "Picornavirus"))   


dat$novel <- as.factor(dat$novel)

# #rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)
# dat.2 <- cbind.data.frame(tip_label = rooted.tree$tip.label)
# dat.plot <- merge(dat.2, dat, by ="tip_label", all.x = T, sort=F)
#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = "none", legend.title = element_blank())
#quartz(p)
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now check that you don't have NAs and blanks throughout the component parts of the name
#you also can edit the original csv file to replace these blanks if you can find the correct info

dat$Isolate  #some are blank, so convert those to NA
dat$Isolate[dat$Isolate==""] <- NA
dat$Accession #all good

#give false numbers to mada samples
#dat$Accession[is.na(dat$Accession)] <- paste(sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",1), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",2), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",3), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",4), sep="_")

dat$Host #no blanks because I manually did them
#dat$Host[dat$Host==""] <- NA
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

#Isolate + accession only
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)] #no

#Isolate + Country only
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] #yes
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
                                                                                                          dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
                                                                                                          #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|", 
                                                                                                          dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|",
                                                                                                          #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                          dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)])

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
#accession Country Host
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
tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = T, sort = F)

names(tree.dat)

tree.dat$tip_label <- tree.dat$new_label
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, bat_Host, Country, Collection_Date, Genus, novel, old_tip_label)

rooted.tree$tip.label <- tree.dat$tip_label

head(tree.dat)
#verify that the old labels match the new
cbind(tree.dat$old_tip_label, rooted.tree$tip.label) #they match

#check out the labels
tree.dat$tip_label#all look good

tree.dat$bat_Host[tree.dat$bat_Host==0] <- "non-bat Host"
tree.dat$bat_Host[tree.dat$bat_Host==1] <- "bat Host"
tree.dat$bat_Host <- as.factor(tree.dat$bat_Host)
shapez = c("bat Host" =  17, "non-bat Host" = 15)
colz2 = c('1' =  "yellow", '0' = "white")


p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=bat_Host), size=3, show.legend = T) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = -.08) +
  scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", family="Helvetica", label.size = 0, alpha=.3, size=2, show.legend=F) +#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=3,y=-5, linesize = .5) +
  theme(legend.position = "left", legend.title = element_blank()) +
  xlim(c(0,4))
p2

##collapse unnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)
#241 megrivirus
#233 mischivirus
#221 cosavirus
#227 cardiovirus
#220 mosavirus
#230 malagasivirus
#212 apthovirus
#208 hunnivirus
#192 salivirus
#197 sicinivirus












##Picornaviridae full and partial rdrp region
setwd(paste0(homewd,"/raxml_trees/picornaviridae_all_rdrp_region"))

#load the tree
tree <-  read.tree("T2.raxml.supportFBP")

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("picornaviridae_rdrp_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #179
length(rooted.tree$tip.label) #179

#check subgroup names
unique(dat$Genus)

colz = c("Cardiovirus" = "cadetblue",    "Enterovirus"  = "cadetblue1",   "Hepatovirus"  = "cadetblue2",   
         "Kobuvirus"   = "cadetblue3" ,   "Parechovirus" = "chartreuse" ,"Erbovirus"  = "chartreuse3" ,    
         "Teschovirus"  = "chartreuse4" ,  "Sapelovirus" = "chocolate" ,   "Tremovirus"  = "chocolate2" ,   
         "Anativirus"   = "chocolate3" ,"Avihepatovirus"  = "chocolate4","Aquamavirus"  = "coral",   
         "Aphthovirus"  = "coral2" ,  "Senecavirus"  = "coral3" ,  "Cosavirus"  = "cyan",
         "Salivirus"   = "cyan3",    "Passerivirus"  = "cyan4", "Oscivirus"   = "darkgoldenrod1" , 
         "Unclassified picornavirus"= "thistle1"  ,   "Mischivirus" = "darkgoldenrod3","Pasivirus"   = "firebrick1"  ,  
         "Gallivirus"   = "firebrick3" ,  "Limnipivirus"  = "firebrick4",  "Hunnivirus"   = "darkolivegreen1",   
         "Dicipivirus" = "darkolivegreen3","Megrivirus"  = "darkorchid1" ,   "Potamipivirus"  = "darkorchid4", 
         "Sakobuvirus"  = "darkseagreen1" ,  "Sicinivirus"  = "darkseagreen3" ,  "Aalivirus"  = "deeppink",
         "Mosavirus"  = "deeppink3" ,    "Rosavirus"  = "deeppink4" ,    "Avisivirus"   = "hotpink",   
         "Orivirus"   = "hotpink3" ,    "Crohivirus" = "hotpink4","Torchivirus" = "indianred1" ,   
         "Bopivirus"  = "indianred3"  ,   "Malagasivirus" = "indianred4" , "Harkavirus"  = "mediumpurple1" ,   
         "Ampivirus"  = "mediumpurple3" ,"Livupivirus" = "mediumpurple4" ,   "Kunsagivirus"  = "orange",  
         "Shanbavirus"  = "orange3" ,  "Rafivirus"   = "seagreen1",  "Coronavirus" ="black",  
         "Poecivirus"  = "seagreen3" ,"Rabovirus"   = "seagreen4",    "Tottorivirus"  = "skyblue1" , 
         "Ailurivirus" = "skyblue3", "Madagascar bat kobuvirus" ="cadetblue3", "Bat picornavirus"="yellow", "Picornavirus"="grey2")

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
                                          "Shanbavirus",  "Rafivirus",  "Coronavirus",  
                                          "Poecivirus", "Rabovirus",    "Tottorivirus", 
                                          "Ailurivirus", "Madagascar bat kobuvirus", "Bat picornavirus", "Picornavirus"))   


dat$novel <- as.factor(dat$novel)

# #rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)
# dat.2 <- cbind.data.frame(tip_label = rooted.tree$tip.label)
# dat.plot <- merge(dat.2, dat, by ="tip_label", all.x = T, sort=F)
#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = "none", legend.title = element_blank())
#quartz(p)
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now check that you don't have NAs and blanks throughout the component parts of the name
#you also can edit the original csv file to replace these blanks if you can find the correct info

dat$Isolate  #some are blank, so convert those to NA
dat$Isolate[dat$Isolate==""] <- NA
dat$Accession #all good

#give false numbers to mada samples
#dat$Accession[is.na(dat$Accession)] <- paste(sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",1), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",2), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",3), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",4), sep="_")

dat$Host #no blanks because I manually did them
#dat$Host[dat$Host==""] <- NA
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

#Isolate + accession only
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)] #no

#Isolate + Country only
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] #yes
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
                                                                                                          dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
                                                                                                          #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|", 
                                                                                                          dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|",
                                                                                                          #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                          dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)])

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
#accession Country Host
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
tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = T, sort = F)

names(tree.dat)

tree.dat$tip_label <- tree.dat$new_label
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, bat_Host, Country, Collection_Date, Genus, novel, old_tip_label)

rooted.tree$tip.label <- tree.dat$tip_label

head(tree.dat)
#verify that the old labels match the new
cbind(tree.dat$old_tip_label, rooted.tree$tip.label) #they match

#check out the labels
tree.dat$tip_label#all look good

tree.dat$bat_Host[tree.dat$bat_Host==0] <- "non-bat Host"
tree.dat$bat_Host[tree.dat$bat_Host==1] <- "bat Host"
tree.dat$bat_Host <- as.factor(tree.dat$bat_Host)
shapez = c("bat Host" =  17, "non-bat Host" = 15)
colz2 = c('1' =  "yellow", '0' = "white")


p3 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=bat_Host), size=3, show.legend = T) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = -.08) +
  scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", family="Helvetica", label.size = 0, alpha=.3, size=2, show.legend=F) +#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=3,y=-5, linesize = .5) +
  theme(legend.position = "left", legend.title = element_blank()) +
  xlim(c(0,4))
p3

##collapse unnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)
#299 salivirus
#302 gallivirus
#309 megrivirus
#322 rosavirus
#246 cardiovirus
#255 cosavirus
#262 apthovirus
#193 enterovirus
#183 hepatovirus










##Picornaviridae full sequences
setwd(paste0(homewd,"/raxml_trees/picornaviridae_full"))

#load the tree
tree <-  read.tree("T2.raxml.supportFBP")

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("picornaviridae_full_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #178
length(rooted.tree$tip.label) #178

#check subgroup names
unique(dat$Genus)

colz = c("Cardiovirus" = "cadetblue",    "Enterovirus"  = "cadetblue1",   "Hepatovirus"  = "cadetblue2",   
         "Kobuvirus"   = "cadetblue3" ,   "Parechovirus" = "chartreuse" ,"Erbovirus"  = "chartreuse3" ,    
         "Teschovirus"  = "chartreuse4" ,  "Sapelovirus" = "chocolate" ,   "Tremovirus"  = "chocolate2" ,   
         "Anativirus"   = "chocolate3" ,"Avihepatovirus"  = "chocolate4","Aquamavirus"  = "coral",   
         "Aphthovirus"  = "coral2" ,  "Senecavirus"  = "coral3" ,  "Cosavirus"  = "cyan",
         "Salivirus"   = "cyan3",    "Passerivirus"  = "cyan4", "Oscivirus"   = "darkgoldenrod1" , 
         "Unclassified picornavirus"= "thistle1"  ,   "Mischivirus" = "darkgoldenrod3","Pasivirus"   = "firebrick1"  ,  
         "Gallivirus"   = "firebrick3" ,  "Limnipivirus"  = "firebrick4",  "Hunnivirus"   = "darkolivegreen1",   
         "Dicipivirus" = "darkolivegreen3","Megrivirus"  = "darkorchid1" ,   "Potamipivirus"  = "darkorchid4", 
         "Sakobuvirus"  = "darkseagreen1" ,  "Sicinivirus"  = "darkseagreen3" ,  "Aalivirus"  = "deeppink",
         "Mosavirus"  = "deeppink3" ,    "Rosavirus"  = "deeppink4" ,    "Avisivirus"   = "hotpink",   
         "Orivirus"   = "hotpink3" ,    "Crohivirus" = "hotpink4","Torchivirus" = "indianred1" ,   
         "Bopivirus"  = "indianred3"  ,   "Malagasivirus" = "indianred4" , "Harkavirus"  = "mediumpurple1" ,   
         "Ampivirus"  = "mediumpurple3" ,"Livupivirus" = "mediumpurple4" ,   "Kunsagivirus"  = "orange",  
         "Shanbavirus"  = "orange3" ,  "Rafivirus"   = "seagreen1",  "Coronavirus" ="black",  
         "Poecivirus"  = "seagreen3" ,"Rabovirus"   = "seagreen4",    "Tottorivirus"  = "skyblue1" , 
         "Ailurivirus" = "skyblue3", "Madagascar bat kobuvirus" ="cadetblue3", "Bat picornavirus"="yellow", "Picornavirus"="grey2")

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
                                          "Shanbavirus",  "Rafivirus",  "Coronavirus",  
                                          "Poecivirus", "Rabovirus",    "Tottorivirus", 
                                          "Ailurivirus", "Madagascar bat kobuvirus", "Bat picornavirus", "Picornavirus"))   


dat$novel <- as.factor(dat$novel)

# #rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)
# dat.2 <- cbind.data.frame(tip_label = rooted.tree$tip.label)
# dat.plot <- merge(dat.2, dat, by ="tip_label", all.x = T, sort=F)
#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = "none", legend.title = element_blank())
#quartz(p)
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now check that you don't have NAs and blanks throughout the component parts of the name
#you also can edit the original csv file to replace these blanks if you can find the correct info

dat$Isolate  #some are blank, so convert those to NA
dat$Isolate[dat$Isolate==""] <- NA
dat$Accession #all good

#give false numbers to mada samples
#dat$Accession[is.na(dat$Accession)] <- paste(sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",1), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",2), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",3), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",4), sep="_")

dat$Host #no blanks because I manually did them
#dat$Host[dat$Host==""] <- NA
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

#Isolate + accession only
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)] #no

#Isolate + Country only
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] #yes
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
                                                                                                          dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
                                                                                                          #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|", 
                                                                                                          dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|",
                                                                                                          #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                          dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)])

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
#accession Country Host
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
tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = T, sort = F)

names(tree.dat)

tree.dat$tip_label <- tree.dat$new_label
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, bat_Host, Country, Collection_Date, Genus, novel, old_tip_label)

rooted.tree$tip.label <- tree.dat$tip_label

head(tree.dat)
#verify that the old labels match the new
cbind(tree.dat$old_tip_label, rooted.tree$tip.label) #they match

#check out the labels
tree.dat$tip_label#all look good

tree.dat$bat_Host[tree.dat$bat_Host==0] <- "non-bat Host"
tree.dat$bat_Host[tree.dat$bat_Host==1] <- "bat Host"
tree.dat$bat_Host <- as.factor(tree.dat$bat_Host)
shapez = c("bat Host" =  17, "non-bat Host" = 15)
colz2 = c('1' =  "yellow", '0' = "white")


p4 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=bat_Host), size=3, show.legend = T) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = -.08) +
  scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", family="Helvetica", label.size = 0, alpha=.3, size=2, show.legend=F) +#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=3,y=-5, linesize = .5) +
  theme(legend.position = "left", legend.title = element_blank()) +
  xlim(c(0,4))
p4

##collapse unnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)
#333 hepatovirus
#324 limnipivirus
#315 parechovirus
#277 salivirus
#270 gallivirus
#296 megrivirus
#301 rosavirus
#190 megrivirus
#349 cosavirus
#200 apthovirus
#237 enterovirus
#236 mosavirus













##Picornaviridae partial sequences
setwd(paste0(homewd,"/raxml_trees/picornaviridae_partial"))

#load the tree
tree <-  read.tree("T2.raxml.supportFBP")

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
#take a quick look in base R
plot(rooted.tree)

#load tree data prepared from elsewhere
dat <- read.csv(("picornaviridae_partial_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check that your metadata matches your tree data
setdiff(rooted.tree$tip.label, dat$tip_label)
#check for duplicates
setdiff(dat$tip_label, rooted.tree$tip.label) #no duplicates
nrow(dat) #172
length(rooted.tree$tip.label) #172

#check subgroup names
unique(dat$Genus)

colz = c("Cardiovirus" = "cadetblue",    "Enterovirus"  = "cadetblue1",   "Hepatovirus"  = "cadetblue2",   
         "Kobuvirus"   = "cadetblue3" ,   "Parechovirus" = "chartreuse" ,"Erbovirus"  = "chartreuse3" ,    
         "Teschovirus"  = "chartreuse4" ,  "Sapelovirus" = "chocolate" ,   "Tremovirus"  = "chocolate2" ,   
         "Anativirus"   = "chocolate3" ,"Avihepatovirus"  = "chocolate4","Aquamavirus"  = "coral",   
         "Aphthovirus"  = "coral2" ,  "Senecavirus"  = "coral3" ,  "Cosavirus"  = "cyan",
         "Salivirus"   = "cyan3",    "Passerivirus"  = "cyan4", "Oscivirus"   = "darkgoldenrod1" , 
         "Unclassified picornavirus"= "thistle1"  ,   "Mischivirus" = "darkgoldenrod3","Pasivirus"   = "firebrick1"  ,  
         "Gallivirus"   = "firebrick3" ,  "Limnipivirus"  = "firebrick4",  "Hunnivirus"   = "darkolivegreen1",   
         "Dicipivirus" = "darkolivegreen3","Megrivirus"  = "darkorchid1" ,   "Potamipivirus"  = "darkorchid4", 
         "Sakobuvirus"  = "darkseagreen1" ,  "Sicinivirus"  = "darkseagreen3" ,  "Aalivirus"  = "deeppink",
         "Mosavirus"  = "deeppink3" ,    "Rosavirus"  = "deeppink4" ,    "Avisivirus"   = "hotpink",   
         "Orivirus"   = "hotpink3" ,    "Crohivirus" = "hotpink4","Torchivirus" = "indianred1" ,   
         "Bopivirus"  = "indianred3"  ,   "Malagasivirus" = "indianred4" , "Harkavirus"  = "mediumpurple1" ,   
         "Ampivirus"  = "mediumpurple3" ,"Livupivirus" = "mediumpurple4" ,   "Kunsagivirus"  = "orange",  
         "Shanbavirus"  = "orange3" ,  "Rafivirus"   = "seagreen1",  "Coronavirus" ="black",  
         "Poecivirus"  = "seagreen3" ,"Rabovirus"   = "seagreen4",    "Tottorivirus"  = "skyblue1" , 
         "Ailurivirus" = "skyblue3", "Madagascar bat kobuvirus" ="cadetblue3", "Bat picornavirus"="yellow", "Picornavirus"="grey2")

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
                                          "Shanbavirus",  "Rafivirus",  "Coronavirus",  
                                          "Poecivirus", "Rabovirus",    "Tottorivirus", 
                                          "Ailurivirus", "Madagascar bat kobuvirus", "Bat picornavirus", "Picornavirus"))   


dat$novel <- as.factor(dat$novel)

# #rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)
# dat.2 <- cbind.data.frame(tip_label = rooted.tree$tip.label)
# dat.plot <- merge(dat.2, dat, by ="tip_label", all.x = T, sort=F)
#take a glance
p <- ggtree(rooted.tree) %<+% dat + geom_tippoint(aes(color=Genus)) +
  geom_tiplab(size=2) + geom_nodelab(size=1) +
  scale_color_manual(values=colz) + theme(legend.position = "none", legend.title = element_blank())
#quartz(p)
p #looks great

#now get new tip labels
dat$old_tip_label <- dat$tip_label
dat$new_label <- NA

#now check that you don't have NAs and blanks throughout the component parts of the name
#you also can edit the original csv file to replace these blanks if you can find the correct info

dat$Isolate  #some are blank, so convert those to NA
dat$Isolate[dat$Isolate==""] <- NA
dat$Accession #all good

#give false numbers to mada samples
#dat$Accession[is.na(dat$Accession)] <- paste(sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",1), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",2), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",3), sapply(strsplit(dat$old_tip_label[is.na(dat$Accession)], "_"), "[",4), sep="_")

dat$Host #no blanks because I manually did them
#dat$Host[dat$Host==""] <- NA
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

#Isolate + accession only
dat$new_label[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)] #no

#Isolate + Country only
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] #yes
dat$new_label[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)] <- paste(dat$Accession[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
                                                                                                          dat$Genus[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|", 
                                                                                                          #dat$Isolate[!is.na(dat$Isolate) & !is.na(dat$Accession) &!is.na(dat$Host) & is.na(dat$Country)], "|", 
                                                                                                          dat$Host[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)], "|",
                                                                                                          #dat$Country[is.na(dat$Isolate) & is.na(dat$Accession) & !is.na(dat$Host) &!is.na(dat$Country)], "|",
                                                                                                          dat$Collection_Date[is.na(dat$Isolate) & !is.na(dat$Accession) & !is.na(dat$Host) &is.na(dat$Country)])

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
#accession Country Host
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
tree.dat <- merge(tree.dat, dat, by = "old_tip_label", all.x = T, sort = F)

names(tree.dat)

tree.dat$tip_label <- tree.dat$new_label
tree.dat <- dplyr::select(tree.dat, tip_label, Isolate, Host, bat_Host, Country, Collection_Date, Genus, novel, old_tip_label)

rooted.tree$tip.label <- tree.dat$tip_label

head(tree.dat)
#verify that the old labels match the new
cbind(tree.dat$old_tip_label, rooted.tree$tip.label) #they match

#check out the labels
tree.dat$tip_label#all look good

tree.dat$bat_Host[tree.dat$bat_Host==0] <- "non-bat Host"
tree.dat$bat_Host[tree.dat$bat_Host==1] <- "bat Host"
tree.dat$bat_Host <- as.factor(tree.dat$bat_Host)
shapez = c("bat Host" =  17, "non-bat Host" = 15)
colz2 = c('1' =  "yellow", '0' = "white")


p5 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=bat_Host), size=3, show.legend = T) +
  geom_nodelab(size=2,nudge_x = -.07, nudge_y = -.08) +
  scale_fill_manual(values=colz) +
  scale_shape_manual(values=shapez) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", family="Helvetica", label.size = 0, alpha=.3, size=2, show.legend=F) +#
  scale_fill_manual(values=colz2) +
  geom_treescale(fontsize=4, x=3,y=-5, linesize = .5) +
  theme(legend.position = "left", legend.title = element_blank()) +
  xlim(c(0,4))
p5

##collapse unnecesary clades
ggtree(rooted.tree) + geom_text(aes(label=node), hjust=-.3)
#318 megrivirus
#314 avisivirus
#251 kobuvirus
#261 salivirus
#245 gallivirus
#268 rosavirus
#273 apthovirus
#179 cosavirus
#340 cardiovirus
#207 enterovirus







