rm(list=ls())

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
#install.packages('gdata')
library(gdata)

homewd= "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus/"


##Picornaviridae all full sequences
setwd(paste0(homewd,"/refseq_raxml_trees/picornaviridae_full"))

#load the tree
tree <-  read.tree("T10.raxml.supportFBP")

rooted.tree <- root(tree, which(tree$tip.label == "NC_030886.1"))
rooted.tree<-drop.tip(rooted.tree, 'NC_045762.1') #drop extra norovirus sequence
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
length(rooted.tree$tip.label) #179

#check subgroup names
unique(dat$Genus)

colz=c("Cardiovirus" = "cadetblue1",    "Enterovirus"  = "cadetblue3",   "Hepatovirus"  = "cadetblue4",  "Kobuvirus"   = "cadetblue2" ,              
       "Parechovirus" = "coral1" ,"Erbovirus"  = "coral2" , "Teschovirus"  = "coral3" ,  "Sapelovirus" = "coral4" , "Tremovirus"  = "goldenrod1" ,
       "Anativirus"   = "darkorange1" ,"Avihepatovirus"  = "darkorange3","Aquamavirus"  = "darkorange4", "Aphthovirus"  = "goldenrod2" ,             
       "Senecavirus"  = "goldenrod3" ,  "Cosavirus"  = "darkslategray1", "Salivirus"   = "darkslategray2",    "Passerivirus"  = "darkslategray3",         
       "Oscivirus"   = "darkslategray4" ,"Unclassified picornavirus"= "indianred1"  ,   "Mischivirus" = "indianred2","Pasivirus"   = "indianred3"  ,
       "Gallivirus"   = "indianred4" ,  "Limnipivirus"  = "lightpink1",  "Hunnivirus"   = "lightpink2",  
       "Dicipivirus" = "lightpink3","Megrivirus"  = "orange1" ,   "Potamipivirus"  = "orange2", "Sakobuvirus"  = "orange1" ,     
       "Sicinivirus"  = "orange3" ,  "Aalivirus"  = "orange4","Mosavirus"  = "orange3" ,    "Rosavirus"  = "pink1" ,          
       "Avisivirus"   = "pink2", "Orivirus"   = "pink3" ,    "Crohivirus" = "plum1","Torchivirus" = "plum2" ,   
       "Bopivirus"  = "plum3"  ,   "Malagasivirus" = "plum4" , "Harkavirus"  = "royalblue2" ,"Ampivirus"  = "royalblue3" ,           
       "Coronavirus" ="black", "Livupivirus" = "royalblue4" ,   "Kunsagivirus"  = "paleturquoise1",  "Shanbavirus"  = "paleturquoise2" ,
       "Rafivirus"   = "paleturquoise3",  "Norovirus"  , "Poecivirus"  = "thistle1" ,"Rabovirus"   = "thistle2",         
       "Tottorivirus"  = "thistle3" , "Ailurivirus" = "yellow1", "Madagascar bat kobuvirus" ="yellow2", "Bat picornavirus"="slategray2",
       "Picornavirus"="slategray3")     



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
                                          "Ailurivirus", "Madagascar bat kobuvirus", "Bat picornavirus", 
                                          "Picornavirus","Coronavirus","Norovirus"))   


dat$novel <- as.factor(dat$novel)

# #rooted.tree.A$node.label <- round(as.numeric(rooted.tree.A$node.label)*100, 0)
# dat.2 <- cbind.data.frame(tip_label = rooted.tree$tip.label)
# dat.plot <- merge(dat.2, dat, by ="tip_label", all.x = T, sort=F)
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

tree.dat$bat_Host[tree.dat$bat_Host==0] <- "Non-bat host"
tree.dat$bat_Host[tree.dat$bat_Host==1] <- "Bat host"
tree.dat$bat_Host <- as.factor(tree.dat$bat_Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")


##uncollapsed tree
p1 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=bat_Host), size=3, show.legend = T) +
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
p2 <- ggtree(rooted.tree) %<+% tree.dat + geom_tippoint(aes(color=Genus, shape=bat_Host), size=3, show.legend = T) +
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
  geom_cladelabel(node = 310, label = "bold(Cardiovirus)",parse=T,hjust='center',offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 297, label = "bold(Cardiovirus)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 282, label = "bold(Megrivirus/Gallivirus)",parse=T,hjust='center', offset.text=0.35, fontsize = 2, color="black") +
  geom_cladelabel(node = 315, label = "bold(Rosavirus)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 272, label = "bold(Apthovirus)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") +
  geom_cladelabel(node = 253, label = "bold(Unclassified)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 242, label = "bold(Avisivirus)",parse=T,hjust='center', offset.text=0.2, fontsize = 2, color="black") +
  geom_cladelabel(node = 238, label = "bold(Parechovirus)",parse=T,hjust='center', offset.text=0.25,fontsize = 2, color="black") +
  geom_cladelabel(node = 247, label = "bold(Limnipivirus)", parse=T,fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 345, label = "bold(Hepatovirus)",parse=T, fontsize = 2,hjust='center', offset.text=0.6, color="black") +
  geom_cladelabel(node = 229, label = "bold(Unclassified)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 188, label = "bold(Enterovirus)",parse=T, fontsize = 2,hjust='center', offset.text=0.2, color="black") +
  geom_cladelabel(node = 184, label = "bold(Enterovirus)",parse=T,hjust='center', offset.text=0.2,fontsize = 2, color="black") 
p2


#collapse the labeled clades
p3<-collapse(p2, 310)+geom_point2(aes(subset=(node==310)), size=3, shape=15, color="cadetblue1", alpha=0.99)
p4<-collapse(p3, 297)+geom_point2(aes(subset=(node==297)), size=3, shape=15, color="cadetblue1")
p5<-collapse(p4, 282)+geom_point2(aes(subset=(node==282)), size=3, shape=15, color="purple")
p6<-collapse(p5, 315)+geom_point2(aes(subset=(node==315)), size=3, shape=15, color="pink1")
p7<-collapse(p6, 272)+geom_point2(aes(subset=(node==272)), size=3, shape=15, color="goldenrod2")
p8<-collapse(p7, 253)+geom_point2(aes(subset=(node==253)), size=3, shape=15, color="indianred1")
p9<-collapse(p8, 242)+geom_point2(aes(subset=(node==242)), size=3, shape=15, color="pink2")
p10<-collapse(p9, 238)+geom_point2(aes(subset=(node==238)), size=3, shape=15, color="coral1")
p11<-collapse(p10, 247)+geom_point2(aes(subset=(node==247)), size=3, shape=15, color="lightpink1")
p12<-collapse(p11, 345)+geom_point2(aes(subset=(node==345)), size=3, shape=15, color="cadetblue4")
p13<-collapse(p12, 229)+geom_point2(aes(subset=(node==229)), size=3, shape=15, color="indianred1")
p14<-collapse(p13, 188)+geom_point2(aes(subset=(node==188)), size=3, shape=15, color="cadetblue3")
p15<-collapse(p14, 184)+geom_point2(aes(subset=(node==184)), size=3, shape=15, color="cadetblue3")
p15

##Clades to collapse
#310 Cardiovirus
#297 cardiovirus
#282 Megrivirus/gallivirus
#315 Rosavirus
#272 Apthovirus
#253 unclassified
#242 avisivirus
#238 Parechovirus
#247 Limnipivirus
#345 hepatovirus
#229 unclassified
#188 enterovirus
#184 enterovirus






