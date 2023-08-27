rm(list=ls())


library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(lubridate)
#library(rBt)
library(treeio)


#make Bayesian timetree from picornaviridae strict molecular clock model

#first, read in the tree

homewd= "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus/"

setwd(paste0(homewd, "/BEAST/refseq_picornaviridae/strict"))

tree <-  read.beast(file = paste0(homewd, "/BEAST/refseq_picornaviridae/strict/picornsaviridae_dtrict_mean_tree"))

treedat <- cbind.data.frame(tip_name = tree@phylo$tip.label)
treedat$beast_name <-treedat$tip_name


treedat$accession_num <- sapply(strsplit(treedat$tip_name, "_"), function(x) x[[1]])
treedat$accession_num[treedat$accession_num=="NC"] <- c("NC_001366","NC_001430","NC_001472","NC_001479","NC_001489","NC_001490",
                                                        "NC_001612","NC_001617","NC_001859","NC_001918","NC_002058","NC_003976",
                                                        "NC_003983","NC_003985","NC_003987","NC_003988","NC_003990","NC_004421",
                                                        "NC_004441","NC_004451","NC_006553","NC_008250","NC_008714","NC_009448",
                                                        "NC_009891","NC_009996","NC_010354","NC_010415","NC_010810","NC_011349",
                                                        "NC_011829","NC_012798","NC_012800","NC_012801","NC_012802","NC_012957",
                                                        "NC_012986","NC_013695","NC_014411","NC_014412","NC_014413","NC_015626",
                                                        "NC_015934","NC_015936","NC_015940","NC_015941","NC_016156", "NC_016403",
                                                        "NC_016769","NC_016964","NC_018226","NC_018400","NC_018506","NC_018668",
                                                        "NC_021178","NC_021201","NC_021220","NC_021482","NC_022332","NC_022802",
                                                        "NC_023162","NC_023422","NC_023858","NC_023861","NC_023984","NC_023985",
                                                        "NC_023987","NC_023988","NC_024070","NC_024073","NC_024120","NC_024765",
                                                        "NC_024766","NC_024767","NC_024768","NC_024769","NC_024770","NC_025114",
                                                        "NC_025432","NC_025474","NC_025675","NC_025890","NC_025961","NC_026249",
                                                        "NC_026315","NC_026316","NC_026470","NC_026921","NC_027054","NC_027214",
                                                        "NC_027818","NC_027918","NC_027919","NC_028363","NC_028364","NC_028365",
                                                        "NC_028366","NC_028380","NC_028964","NC_028970","NC_028981","NC_029854",
                                                        "NC_029905","NC_030454","NC_030843","NC_031105","NC_031106","NC_032126",
                                                        "NC_033695","NC_033793","NC_033818","NC_033819","NC_033820","NC_034206",
                                                        "NC_034245","NC_034267","NC_034381","NC_034385","NC_034453","NC_034617",
                                                        "NC_034971","NC_035110","NC_035198","NC_035779","NC_036588","NC_037654",
                                                        "NC_038303","NC_038304","NC_038305","NC_038306","NC_038307","NC_038308",
                                                        "NC_038309","NC_038310","NC_038311","NC_038312","NC_038313","NC_038314",
                                                        "NC_038315","NC_038316","NC_038317","NC_038318","NC_038319","NC_038878",
                                                        "NC_038880","NC_038957","NC_038961","NC_038989","NC_039004","NC_039209",
                                                        "NC_039210","NC_039212","NC_039235","NC_040605","NC_040611","NC_040642",
                                                        "NC_040673","NC_040684","NC_043072","NC_043544","NC_055108","NC_055156",
                                                        "NC_055159","NC_055160","NC_055161")
treedat$accession_num[treedat$accession_num=="F"] <- c("OQ818316","OQ818317","OQ818318","OQ818320","OQ818321","OQ818322",
                                                       "OQ818323","OQ818324","OQ818325","OQ818328","OQ818329","OQ818337")
treedat$accession_num[treedat$accession_num=="F"] <- c("OP287812")


#and load data of corresponding tree

dat <- read.csv(file = "picornaviridae_beast_metadata.csv", header = T, stringsAsFactors = F)
dat$Collection_Date <- as.Date(dat$Collection_Date)


ccolz = c("Cardiovirus" = "cadetblue1",    "Enterovirus"  = "cadetblue2",   "Hepatovirus"  = "cadetblue3",   
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


mrsd.dat <- max(dat$Collection_Date)
p1 <- ggtree(tree, mrsd=mrsd.dat)  + theme_tree2()  +geom_nodelab()
p1

tree.dat <- p1$data
node.sub <- dplyr::select(tree.dat, node, x)
names(node.sub) <-  c("node", "nodetime")

#and 
head(dat)

dat.plot <- merge(treedat, dat, by="accession_num", all.x = T, sort=F)

head(dat.plot)
dat.plot$new_label = NA
dat.plot$new_label[!is.na(dat.plot$Species)] <- paste(dat.plot$accession_num[!is.na(dat.plot$Species)], " | ", 
                                                     dat.plot$Species[!is.na(dat.plot$Species)], " | ", 
                                                     dat.plot$source[!is.na(dat.plot$Species)], " | ",
                                                     dat.plot$country[!is.na(dat.plot$Species)], " | ",
                                                     dat.plot$collection_year[!is.na(dat.plot$Species)])

dat.plot$new_label[is.na(dat.plot$Species)] <- paste(dat.plot$accession_num[is.na(dat.plot$Species)], " | ", 
                                                    dat.plot$source[is.na(dat.plot$Species)], " | ",
                                                    dat.plot$country[is.na(dat.plot$Species)], " | ",
                                                    dat.plot$collection_year[is.na(dat.plot$Species)])


tree@phylo$tip.label <- dat.plot$new_label

dat.sub <- dplyr::select(dat.plot, new_label, Genus, Family, Collection_Date,Host, Country, source, Collection_Year, novel, Isolate, Species)
head(dat.sub)

dat.sub$Host[dat.sub$Host==0] <- "Non-bat host"
dat.sub$Host[dat.sub$Host==1] <- "Bat host"
dat.sub$Host <- as.factor(dat.sub$Host)
shapez = c("Bat host" =  17, "Non-bat host" = 19)
colz2 = c('1' =  "yellow", '0' = "white")

#plot tree
p1 <-ggtree(tree, mrsd=mrsd.dat, size=.8) %<+% dat.sub+
  geom_tippoint(aes(color=Genus, shape=Host), size=3,stroke=0,show.legend = T) +  
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz)+
  scale_shape_manual(values=shapez) +
  new_scale_fill() +
  theme_tree2() +
  coord_cartesian(clip = "off") + 
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=3, stroke=.1)+
  scale_fill_continuous(low="yellow", high="red") +
  guides(fill_continuous = guide_legend(order = 2),col = guide_legend(order = 1)) +
  theme(legend.position = "left", 
        legend.direction = "vertical",
        legend.text = element_text(size=8), 
        legend.key.size = unit(0.3, "cm"))+
  scale_x_continuous(breaks=c(1400, 1600, 1800, 2000),
                     labels=c(623, 423, 223, 23)) +
  xlab("Years to MRCA")+ggnewscale::new_scale_fill() + 
  geom_tiplab(aes(fill = novel, show.legend=F), geom = "label", Genus="Helvetica", label.size = 0, label.padding = unit(0, "lines"), alpha=.4, size=3,  nudge_x=0.1) +
  scale_fill_manual(values=colz2) +
  guides(fill="none")
p1




