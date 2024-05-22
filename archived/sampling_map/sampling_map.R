#This is modified from the Mada-Bat-Cov-main GitHub repo for figure 1. This script will make a figure of 
#a map of Madagascar that shows breakdown of sample metadata with percent positives and negatives

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

# To run this script, change the "mainwd" to wherever this folder
# ("Mada_bat_picorna") is stored on your computer
# Also, make sure to download/clone the "Mada-GIS" folder to 
# your home computer. I recommend putting it in the same parent 
# directory as "Mada_bat_picorna". My files are on the desktop so make
# change your directory to where the files are. 

# For example, my two folders are stored at:

# "/Users/gwenddolenkettenburg/Desktop/Mada_bat_picorna/     ...AND
# "/Users/gwenddolenkettenburg/Desktop/Mada-GIS-main/

# I keep all my github repos under "R_repositories"

#####################################################################
#####################################################################
# Set wd to data on this computer. Also ID homewd, assuming that 
# Mada-GIS is cloned to the same series of sub-folders
homewd = "/Users/gwenddolenkettenburg/Desktop/mada-bat-picornavirus" 
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
#print(p2)

  ggsave(file = "general_map_compass_sampling_locations.pdf",
         plot = p2,
          units="mm",  
          width=40, 
          height=60, 
          scale=3, 
          dpi=300)
# # 
coordinate$label <- coordinate$species
coordinate$label[coordinate$label=="Pteropus rufus"] <- "Pteropus\nrufus"
coordinate$label[coordinate$label=="Rousettus madagascariensis"] <- "Rousettus\nmadagascariensis"
coordinate$label[coordinate$label=="Eidolon dupreanum"] <- "Eidolon\ndupreanum"


#load GPS point and label
p2b<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="black",size=1,data=coordinate)+
  geom_text(data= coordinate,                       #### Labeling
            aes(x=longitude_e, y=latitude_s, label=label),
            fontface="italic",
            color = "#1B262C", size=3,
            nudge_x = c(4.5,2,5.5),
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
# # # 
     ggsave(file = "general_map_with_names_locations.pdf",
            plot = p2b,
            units="mm",  
            width=40, 
            height=60, 
            scale=3, 
            dpi=300)
# # # 


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
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.8,.8),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.text = element_text(size = 7.5)) +
  scale_fill_manual(values=colz) +
  geom_scatterpie_legend(log10(c(10,100)/1.2),
                         x=54, y=-23.5, 
                         n=2,
                         labeller = function(x) paste(10^(x)*1,"indiv"))

p4

