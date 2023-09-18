###This is a bunch of separate scripts pasted together to make figure 1

rm(list=ls())

#packages
library(sf)
library(mapplots)
library(scatterpie)
library(maptools)
library(plyr) 
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(ggspatial)
library(ggrepel)
library(hrbrthemes)
library(reshape)
library(paletteer)
library(ggh4x)
library(cowplot)
library(ggplotify)

##Part A: map

# Set wd to data on this computer. Also ID homewd, assuming that 
# Mada-GIS is cloned to the same series of sub-folders
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 
#should be wherever "Mada_bat_picorna" is stored on your home computer
#basewd = paste(strsplit(homewd, "/")[[1]][1:6], collapse = "/")
#Cara had the above command as basewd but it did not work for me and kept going into the Mada_bat_picorna folder
#So using Sophia's astrovirus map file I changed it to the basewd command below
basewd = "/Users/gwenddolenkettenburg/Desktop/Mada-GIS"
#mapwd = paste0(basewd, "/", "Mada_GIS")
setwd(paste0(homewd, "/", "sampling_map/"))

#import madagascar shapfile
name<- paste0(basewd, "/", "MDG-3/MDG_adm3.shp")
otl_file <- paste(name, sep="") 
orotl_shp <- st_read(otl_file)
#View(orotl_shp)  # Open attribute table
class(orotl_shp)

###import and configuration

p1<-ggplot() +  
  geom_sf(color = "coral1", fill = "coral1",data = orotl_shp)+
  coord_sf(xlim = c(42, 60), ylim = c(-26, -11.5), expand = FALSE)+
  theme_bw()+
  theme(plot.margin = unit(c(-1,.5,-1.5,.1),"cm"))+
  xlab("Longitude") + ylab("Latitude") 
#print(p1)


#import picorna data
dat <- read.csv(file = paste0(homewd,"/metadata/demo_picornavirales_fecesurine_meta.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

#clean class
unique(dat$bat_sex)

dat$bat_sex <- dat$bat_sex
dat$bat_sex[dat$bat_sex=="female"] <- "F"
dat$bat_sex[dat$bat_sex=="male"] <- "M"

#clean class
unique(dat$bat_age_class)

#and rank by rough age
unique(dat$young_of_year)
dat$bat_age_class <- dat$bat_age_class
dat$bat_age_class[dat$bat_age_class=="P" | dat$bat_age_class=="L"] <- "A"
dat$bat_age_class[dat$bat_age_class=="NL" | dat$young_of_year=="no"] <- "A"
dat$bat_age_class[dat$young_of_year=="yes"] <- "J"


unique(dat$any_picorna)


# now subset the data to just include the columns of interest
dat <- dplyr::select(dat,roost_site,latitude_s, longitude_e,
                     collection_date,
                     species, sampleid, any_picorna)

head(dat)
unique(dat$roost_site)

#get sites
coordinate <- ddply(dat, .(roost_site), summarise, latitude_s=unique(latitude_s), longitude_e=unique(longitude_e))
coordinate <-subset(coordinate, roost_site=="Ambakoana" | roost_site=="AngavoKely" | roost_site=="Maromizaha")
coordinate$species <- c("Pteropus rufus", "Eidolon dupreanum","Rousettus madagascariensis")
head(coordinate)

#plot sites on map
p2<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="black",size=1,data=dat)+
  annotation_scale(location = "bl", width_hint = 0.05) +    # scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.02, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)

coordinate$label <- coordinate$species
coordinate$label[coordinate$label=="Pteropus rufus"] <- "Pteropus\nrufus"
coordinate$label[coordinate$label=="Rousettus madagascariensis"] <- "Rousettus\nmadagascariensis"
coordinate$label[coordinate$label=="Eidolon dupreanum"] <- "Eidolon\ndupreanum"


#load GPS point and label
p2b<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="black",size=1,data=coordinate)+
  geom_text(data= coordinate,                       #### Labeling
            aes(x=longitude_e, y=latitude_s, label=label),
            fontface="italic",
            color = "#1B262C", size=4,
            nudge_x = c(4.8,2.6,6),
            nudge_y = c(3,-4,-0.3),
            check_overlap = T)+
  annotation_scale(location = "bl", width_hint = 0.05) +    #scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.03, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)+
  geom_text_repel(segment.colour="black")+
  theme_bw() +theme(panel.grid = element_blank(), 
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(-1,.5,-1.5,.1),"cm"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.26,.90),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.background = element_rect(color="gray",size = .1),
                    legend.text = element_text(size = 9,face = "italic"))
p2b


#plot one site per species
dat$roost_site[dat$species=="Pteropus rufus"] <- "Ambakoana"
dat$roost_site[dat$species=="Eidolon dupreanum"] <- "AngavoKely"
dat$roost_site[dat$species=="Rousettus madagascariensis"] <- "Maromizaha"

dat$longitude_e[dat$roost_site=="Ambakoana"] <- coordinate$longitude_e[coordinate$roost_site=="Ambakoana"]
dat$longitude_e[dat$roost_site=="AngavoKely"] <- coordinate$longitude_e[coordinate$roost_site=="AngavoKely"]
dat$longitude_e[dat$roost_site=="Maromizaha"] <- coordinate$longitude_e[coordinate$roost_site=="Maromizaha"]


dat$latitude_s[dat$roost_site=="Ambakoana"] <- coordinate$latitude_s[coordinate$roost_site=="Ambakoana"]
dat$latitude_s[dat$roost_site=="AngavoKely"] <- coordinate$latitude_s[coordinate$roost_site=="AngavoKely"]
dat$latitude_s[dat$roost_site=="Maromizaha"] <- coordinate$latitude_s[coordinate$roost_site=="Maromizaha"]


###by total positives is anything with dat
dat$plot_class<-NA
dat$plot_class[dat$any_picorna==1]<-"Positive"
dat$plot_class[dat$any_picorna==0]<-"Negative"


##by total positives
pies <- ddply(dat, .(species, roost_site, latitude_s, longitude_e, plot_class), summarise, value=length(sampleid))


#tot_sum = ddply(pies,.(species, bat_sex), summarise,N=sum(value))
tot_sum = ddply(pies,.(species), summarise,N=sum(value))

#pies <- merge(pies, tot_sum, by=c("species", "bat_sex"), all.x=T)
pies <- merge(pies, tot_sum, by=c("species"), all.x=T)

#pies$plot_class <- factor(pies$plot_class, levels=c( "male -", "female -", "male +", "female +"))
pies$plot_class <- factor(pies$plot_class, levels=c( "Positive", "Negative"))

#now split into two pies
piesP = subset(pies, plot_class=="Positive")
piesN = subset(pies, plot_class=="Negative")

###Get the pie data in the right format###
colz = c('Positive' ="darkgoldenrod2", 'Negative' ="lightskyblue2")



p3<-ggplot() +
  geom_scatterpie(aes(x=longitude_e, y=latitude_s, r=(N/1000)),
                  data = pies, cols="plot_class", long_format=TRUE) +
  scale_fill_manual(values=colz)

p3

# # copy of latitude (x.) and longitude (y.)
pies$x2 <- pies$longitude_e
pies$y2 <- pies$latitude_s


# #manually move the pie chart in case there is an overlap (change x and y)

#from the labels
# nudge_x = c(4.5 PR,2 ED,5.5 RM),
# nudge_y = c(3 PR,-4 ED,-0.3 RM),

pies$x2[pies$species== "Pteropus rufus"] <- pies$longitude_e[pies$species== "Pteropus rufus"] + 2
pies$y2[pies$species== "Pteropus rufus"] <- pies$latitude_s[pies$species== "Pteropus rufus"] + 3

pies$x2[pies$species== "Eidolon dupreanum"] <- pies$longitude_e[pies$species== "Eidolon dupreanum"] -1
pies$y2[pies$species== "Eidolon dupreanum"] <- pies$latitude_s[pies$species== "Eidolon dupreanum"] - 4

pies$x2[pies$species== "Rousettus madagascariensis"] <- pies$longitude_e[pies$species== "Rousettus madagascariensis"] + 2 
pies$y2[pies$species== "Rousettus madagascariensis"] <- pies$latitude_s[pies$species== "Rousettus madagascariensis"] - 0.3

head(pies)



#This is part A of figure 1
p4 <- p2b+
  annotate("segment", x=pies$longitude_e, xend=pies$x2,y=pies$latitude_s,yend=pies$y2,size=.7)+ # put the lines
  geom_scatterpie(aes(x=x2, y=y2, r=(log10(N)/1.2)), 
                  data = pies, cols="plot_class", long_format=TRUE) +
  theme_bw() +theme(panel.grid = element_blank(),
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(-1,.5,-1.5,.1),"cm"),
                    axis.title.x = element_text(color="black", size=10),
                    axis.title.y = element_text(color="black", size=10),
                    legend.position=c(.8,.8),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.text = element_text(size = 9)) +
  scale_fill_manual(values=colz) +
  geom_scatterpie_legend(log10(c(10,100)/1.2),
                         x=54, y=-23.5, 
                         n=2,
                         labeller = function(x) paste(10^(x)*1,"indiv"))

p4










##Part B: summary heat map of types of picornavirales
dat <- read.csv(file = paste0(homewd,"/metadata/summary_demo_heatmap.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

dat$virus<-factor(dat$virus, levels=c("Bat picornavirus", "Hepatovirus", "Kobuvirus", "Kunsagivirus", 
                                      "Mischivirus", "Sapelovirus", "Teschovirus", "Sapovirus", "Aparavirus",
                                      "Cripavirus", "Cheravirus", "Nepovirus", "Felisavirus", "Picorna-like virus", "Tetnovirus"))
dat$family<-factor(dat$family, levels=c("Picornaviridae", "Caliciviridae", "Dicistroviridae", "Secoviridae", "Unclassified"))


ED<-subset(dat, species=="E. dupreanum")
PR<-subset(dat, species=="P. rufus")
RM<-subset(dat, species=="R. madagascariensis")


picorna<-subset(dat,picorna=="1")
other<-subset(dat,picorna=="0")


#facet by viral family
p9 <- ggplot(dat, aes(x=species, y=virus)) +
  geom_tile(aes(fill = cut(num_genome,breaks=0:3, labels=1:3))) +
  scale_fill_manual(values=c("lightblue1","skyblue", "royalblue1")) +
  #scale_fill_viridis_c(option = "G", direction = -1) +
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  facet_nested(family+type~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  geom_text(aes(label=num_genome))+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Count",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  scale_x_discrete(position="top")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(face="italic"),
        legend.position = "right")
p9





p6 <- ggplot(ED, aes(x=type, y=virus)) +
  geom_tile(aes(fill = cut(num_genome,breaks=0:3, labels=1:3))) +
  scale_fill_manual(values=c("lightblue1","royalblue1")) +
  #scale_fill_viridis_c(option = "G", direction = -1) +
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  facet_nested(family~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  geom_text(aes(label=num_genome))+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Count",
    title="Eidolon dupreanum")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(0, 0, 0, 2),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=9,face="italic"),
        axis.text.x = element_text(size=9),
        legend.position = "none")
p6



p7 <- ggplot(PR, aes(x=type, y=virus)) +
  geom_tile(aes(fill = cut(num_genome,breaks=0:3, labels=1:3))) +
  scale_fill_manual(values=c("lightblue1","royalblue1")) +
  #scale_fill_viridis_c(option = "G", direction = -1) +
  facet_nested(family~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  geom_text(aes(label=num_genome))+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "Genome type",
    y= "",
    fill="Count",
    title="Pteropus rufus")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=9,face="italic"),
        axis.text.x = element_text(size=9),
        legend.position="none")
p7


p8 <- ggplot(RM, aes(x=type, y=virus)) +
  geom_tile(aes(fill = cut(num_genome,breaks=0:3, labels=1:3))) +
  scale_fill_manual(values=c("lightblue1","lightskyblue", "royalblue1")) +
  #scale_fill_viridis_c(option = "G", direction = -1) +
  facet_nested(family~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  geom_text(aes(label=num_genome))+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Count",
    title="Rousettus madagascariensis")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=9,face="italic"),
        axis.text.x = element_text(size=9),
        legend.position=c(0.85,0.9),
        legend.margin = margin(c(0,0,0,0)))
p8

###Put stuff together

##Version 1

final<- plot_grid(p6,p7,p8,labels=c("","",""),rel_widths = c(1, 1,1), rel_heights = c(1, 1,1),ncol=3, align="hv", axis="b")
final

final2<- plot_grid(p4,final,labels=c("A","B"),rel_widths = c(1, 1.8), rel_heights = c(3, 1),ncol=2, align="hv", axis="b", label_size = 23)
final2







#Supplementary info if needed

##Try making a heatmap of the sampling data like age and percent of each bat species positive
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 
dat <- read.csv(file = paste0(homewd,"/metadata/demo_picornavirales_fecesurine_meta_2018_2019.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

#clean class
unique(dat$bat_sex)

dat$bat_sex <- dat$bat_sex
dat$bat_sex[dat$bat_sex=="female"] <- "Female"
dat$bat_sex[dat$bat_sex=="male"] <- "Male"

dat$month[dat$month=="1"]<-"Jan"
dat$month[dat$month=="2"]<-"Feb"
dat$month[dat$month=="3"]<-"March"
dat$month[dat$month=="4"]<-"April"
dat$month[dat$month=="5"]<-"May"
dat$month[dat$month=="6"]<-"June"
dat$month[dat$month=="7"]<-"July"
dat$month[dat$month=="8"]<-"August"
dat$month[dat$month=="9"]<-"Sept"
dat$month[dat$month=="10"]<-"Oct"
dat$month[dat$month=="11"]<-"Nov"
dat$month[dat$month=="12"]<-"Dec"

#clean class
unique(dat$bat_age_class)

#and rank by rough age
unique(dat$young_of_year)
dat$bat_age_class <- dat$bat_age_class
dat$bat_age_class[dat$bat_age_class=="P" | dat$bat_age_class=="L"] <- "A"
dat$bat_age_class[dat$bat_age_class=="NL" | dat$young_of_year=="no"] <- "A"
dat$bat_age_class[dat$young_of_year=="yes"] <- "J"


unique(dat$any_picorna)

#select columns of interest
dat <- dplyr::select(dat,roost_site,bat_age_class,bat_sex, collection_date,
                     species, sampleid, any_picorna, month, year, day_of_year)
dem.dat<-dplyr::select(dat,roost_site,bat_age_class,bat_sex,
                       species, sampleid, any_picorna)

head(dat)
unique(dat$roost_site)

#get into date form
dat$collection_date <- as.Date(dat$collection_date,format = "%m/%d/%y")

#change picorna to numeric
dat$any_picorna <- as.numeric(dat$any_picorna)
names(dat)[names(dat)=="bat_species"] <- "species"

dem.dat$any_picorna <- as.numeric(dem.dat$any_picorna)
names(dem.dat)[names(dem.dat)=="bat_species"] <- "species"

#and make sure it is only 1 of the same sample type per each date
dat.list <- dlply(dat, .(sampleid))
dem.dat.list <- dlply(dem.dat, .(sampleid))

#summarize into prevalence by species and epiwk
dat.sum <- ddply(dat, .(species, month, bat_age_class, bat_sex), summarise, N=length(any_picorna), pos=sum(any_picorna))
dem.dat.sum <- ddply(dem.dat, .(species, bat_age_class, bat_sex), summarise, N=length(any_picorna), pos=sum(any_picorna))

#get negatives and prevalence
dat.sum$neg= dat.sum$N-dat.sum$pos
dat.sum$prevalence <- dat.sum$pos/dat.sum$N

dem.dat.sum$neg= dem.dat.sum$N-dem.dat.sum$pos
dem.dat.sum$prevalence <- dem.dat.sum$pos/dem.dat.sum$N

#and confidence intervals on the prevalence
CIs <- mapply(FUN=prop.test, x=as.list(dat.sum$pos), n=as.list(dat.sum$N), MoreArgs = list(alternative = "two.sided", conf.level = .95, correct=F), SIMPLIFY = F)
CIs <- mapply(FUN=prop.test, x=as.list(dem.dat.sum$pos), n=as.list(dem.dat.sum$N), MoreArgs = list(alternative = "two.sided", conf.level = .95, correct=F), SIMPLIFY = F)


#and extract the upper and lower CIs
get.CIs <- function(df){
  lci = df$conf.int[1]
  uci = df$conf.int[2]
  out.df<- cbind.data.frame(lci=lci, uci=uci)
  return(out.df)
}

CIs <- lapply(CIs, get.CIs)

dat.sum$lci <- c(unlist(sapply(CIs, '[',1)))
dat.sum$uci <- c(unlist(sapply(CIs, '[',2)))

#simplify= the name of "host_genus_species" 
names(dat.sum)[names(dat.sum)=="host_genus_species"] <- "species"
names(dem.dat.sum)[names(dem.dat.sum)=="host_genus_species"] <- "species"


#make the data into factors
dat.sum$species<-factor(dat.sum$species, levels=c("Eidolon dupreanum", "Pteropus rufus", "Rousettus madagascariensis"))
dat.sum$bat_sex<-factor(dat.sum$bat_sex, levels=c("Male", "Female"))
dat.sum$bat_age_class<-factor(dat.sum$bat_age_class, levels=c("J", "A"))
dat.sum$month<-factor(dat.sum$month, levels=c("Jan", "Feb", "March", "April", "May", "June", "July", "August", "Sept",
                                              "Oct", "Nov", "Dec"))

dem.dat.sum$species<-factor(dem.dat.sum$species, levels=c("Eidolon dupreanum", "Pteropus rufus", "Rousettus madagascariensis"))
dem.dat.sum$bat_sex<-factor(dem.dat.sum$bat_sex, levels=c("Male", "Female"))
dem.dat.sum$bat_age_class<-factor(dem.dat.sum$bat_age_class, levels=c("J", "A"))

#dat.sum$any_picorna<-factor(dat.sum$any_picorna, levels=c("0", "1"))

#subset some data
ED1<-subset(dat.sum, species=="Eidolon dupreanum")
PR1<-subset(dat.sum, species=="Pteropus rufus")
RM1<-subset(dat.sum, species=="Rousettus madagascariensis")

ED2<-subset(dem.dat.sum, species=="Eidolon dupreanum")
PR2<-subset(dem.dat.sum, species=="Pteropus rufus")
RM2<-subset(dem.dat.sum, species=="Rousettus madagascariensis")



#facet by month 
p10 <- ggplot(dat.sum, aes(x=species, y=bat_age_class, fill=pos)) +
  geom_tile() +
  #scale_fill_manual(values=c("lightblue1""royalblue1")) +
  scale_fill_viridis_c(option = "F", direction = -1) +
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  facet_nested(month+bat_sex~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Bats \n positive",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(face="italic"),
        legend.position = "right")
p10


#just demographic info
p11 <- ggplot(dem.dat.sum, aes(x=bat_sex, y=bat_age_class, fill=pos)) +
  geom_tile() +
  #scale_fill_manual(values=c("lightblue1""royalblue1")) +
  scale_fill_viridis_c(option = "F", direction = -1) +
  facet_wrap(~species, ncol=3,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  # facet_nested(month+bat_sex~., scales="free", space="free",
  #              switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Bats \n positive",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="right")+
  theme(plot.margin = margin(2, 0, 0, 0),
        legend.margin = margin(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(face="italic"),
        legend.position = "bottom")
p11


###Put stuff together

##Version 2

final<- plot_grid(p6,p7,p8,labels=c("C","",""), rel_widths = c(1, 1,1),hjust=1.5, rel_heights = c(1, 1,1),
                  ncol=3, align="hv", axis="l", label_size=23)
final

final3<-plot_grid(p4,p11, labels = c("A", "B"), rel_widths = c(1,0.8),hjust=0.5, rel_heights = c(1,0.7), ncol=1, align="hv", axis="l", label_size=23)
final3

final4<- plot_grid(NULL,final3,NULL,final,labels=c("","", "",""),rel_widths = c(0.1,1,0.1, 2), rel_heights = c(1,1,1, 0.8),ncol=4, align="hv", axis="b", label_size = 23)
final4


#export landscape PDF inches 20x9 inches