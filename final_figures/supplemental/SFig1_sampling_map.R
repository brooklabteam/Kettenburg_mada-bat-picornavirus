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
library(reshape)
library(paletteer)
library(hrbrthemes)
library(reshape)

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
homewd = "/Users/gwenddolenkettenburg/Desktop/developer/Kettenburg_mada-bat-picornavirus" 
#should be wherever "Mada_bat_picorna" is stored on your home computer
#basewd = paste(strsplit(homewd, "/")[[1]][1:6], collapse = "/")
#Cara had the above command as basewd but it did not work for me and kept going into the Mada_bat_picorna folder
#So using Sophia's astrovirus map file I changed it to the basewd command below
basewd = "/Users/gwenddolenkettenburg/Desktop/developer/Mada-GIS"
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
  geom_sf(color = "sienna4", fill = "sienna1",data = orotl_shp)+
  coord_sf(xlim = c(42, 56), ylim = c(-26, -11.5), expand = FALSE)+
  theme_bw()+
  theme(plot.margin = unit(c(-1,.5,-1.5,.1),"cm"))+
  xlab("Longitude") + ylab("Latitude") 
#print(p1)

  
#import picorna data
dat <- read.csv(file = paste0(homewd,"/metadata/demo_sampleset_fecesurine_map.csv"), header = T, stringsAsFactors = F)
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
                       processing_date,
                       bat_species, sample_id, any_picorna, any_calici, any_pos)

head(dat)
unique(dat$roost_site)

#get sites


#Ankarana does not have any positives so I combined both the E dupreanum and R. madagascariensis just for the purpose of only haveing one pie
coordinate <- ddply(dat, .(roost_site), summarise, latitude_s=unique(latitude_s), longitude_e=unique(longitude_e))
#coordinate <-subset(coordinate, roost_site=="Ambakoana" | roost_site=="Angavokely" | roost_site=="Ankarana_ED" | roost_site=="Ankarana_RM"| roost_site=="Maromizaha")
coordinate <-subset(coordinate, roost_site=="Ambakoana" | roost_site=="Angavokely/Angavobe" | roost_site=="Ankarana" | roost_site=="Maromizaha")
#coordinate$bat_species <- c("Pteropus rufus", "Eidolon dupreanum,"Eidolon dupreanum","Rousettus madagascariensis","Rousettus madagascariensis")
coordinate$bat_species <- c("Pteropus rufus", "Eidolon dupreanum","Eidolon dupreanum/Rousettus madagascariensis","Rousettus madagascariensis")
head(coordinate)

#plot sites on map
p2<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="red",size=1,data=dat)+
  annotation_scale(location = "bl", width_hint = 0.05) +    # scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.02, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)
#print(p2)

coordinate$label <- coordinate$bat_species
coordinate$label[coordinate$label=="Pteropus rufus"] <- "Pteropus\nrufus (N=146)"
coordinate$label[coordinate$label=="Rousettus madagascariensis"] <- "Rousettus\nmadagascariensis \n(N=225)"
coordinate$label[coordinate$label=="Eidolon dupreanum"] <- "Eidolon\ndupreanum \n(N=281)"
coordinate$label[coordinate$label=="Eidolon dupreanum/Rousettus madagascariensis"] <- "N=92 Eidolon dupreanum &\nN=59 Rousettus madagascariensis \n(0% picorna & calici prevalence)"

#order is Ambakoana, Angavokely, Ankarana,  Maromizaha

#Pies are 
# nudge_x = c(6 PR,-1 ED,2 RM,3 ED_RM),
# nudge_y = c(1 PR,-4.4 ED,2 RM, -1 ED_RM),

#nudge_x = c(PR 6,ED 3.5,ED_RM 3.5,RM 5),
#nudge_y = c(PR 1.5,ED -5,ED_RM 0,RM -1.5),

#load GPS point and label
p2b<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="black",size=2,data=coordinate)+
  geom_text(data= coordinate,                       #### Labeling
            aes(x=longitude_e, y=latitude_s, label=label),
            fontface="italic",
            color = "#1B262C", size=3,
            nudge_x = c(6,2.8,3.5,5.4),
            nudge_y = c(1.5,-5,0,-1.5),
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

###by total positives is anything with dat
dat$plot_class_pos<-NA
dat$plot_class_pos[dat$any_pos==1]<-"Picorna positive"
dat$plot_class_pos[dat$any_pos==0]<-"Negative"
dat$plot_class_pos[dat$any_pos==2]<-"Calici positive"


##by total positives

#get rid of ankarana site because it has no positives
drop.dat <- droplevels(dat[!dat$roost_site == 'Ankarana',])

pies_pos <- ddply(drop.dat, .(bat_species, roost_site, latitude_s, longitude_e, plot_class_pos), summarise, value=length(sample_id))


#tot_sum = ddply(pies,.(bat_species, bat_sex), summarise,N=sum(value))
tot_sum_pos = ddply(pies_pos,.(bat_species), summarise,N=sum(value))

#pies <- merge(pies, tot_sum, by=c("bat_species", "bat_sex"), all.x=T)
pies_pos <- merge(pies_pos, tot_sum_pos, by=c("bat_species"), all.x=T)

#pies$plot_class <- factor(pies$plot_class, levels=c( "male -", "female -", "male +", "female +"))
pies_pos$plot_class_pos <- factor(pies_pos$plot_class_pos, levels=c("Picorna positive","Calici positive", "Negative"))

#now split into two or three pies
piesPP_pos = subset(pies_pos, plot_class_pos=="Picorna positive")
piesPC_pos = subset(pies_pos, plot_class_pos=="Calici positive")
piesN_pos = subset(pies_pos, plot_class_pos=="Negative")

###Get the pie data in the right format###
colz = c('Positive' ="goldenrod1", 'Negative' ="lightcyan2", "Picorna positive"="goldenrod1", "Calici positive"="purple1")

p3_pos<-ggplot() +
  geom_scatterpie(aes(x=longitude_e, y=latitude_s, r=(N/1000)),
                  data = pies_pos, cols="plot_class_pos", long_format=TRUE) +
  scale_fill_manual(values=colz)

p3_pos
 
# # copy of latitude (x.) and longitude (y.)
pies_pos$x2 <- pies_pos$longitude_e
pies_pos$y2 <- pies_pos$latitude_s

coordinate$x2<-coordinate$longitude_e
coordinate$y2<-coordinate$latitude_s


# #manually move the pie chart in case there is an overlap (change x and y)

#from the labels
#nudge_x = c(PR 6,ED 3.5,ED_RM 3.5,RM 5),
#nudge_y = c(PR 1.5,ED -5,ED_RM 0,RM -1.5),

#pos
pies_pos$x2[pies_pos$bat_species== "Pteropus rufus"] <- pies_pos$longitude_e[pies_pos$bat_species== "Pteropus rufus"] + 3.3
pies_pos$y2[pies_pos$bat_species== "Pteropus rufus"] <- pies_pos$latitude_s[pies_pos$bat_species== "Pteropus rufus"] + 1.2

pies_pos$x2[pies_pos$bat_species== "Eidolon dupreanum"] <- pies_pos$longitude_e[pies_pos$bat_species== "Eidolon dupreanum"] -1
pies_pos$y2[pies_pos$bat_species== "Eidolon dupreanum"] <- pies_pos$latitude_s[pies_pos$bat_species== "Eidolon dupreanum"] - 4.5

pies_pos$x2[pies_pos$bat_species== "Rousettus madagascariensis"] <- pies_pos$longitude_e[pies_pos$bat_species== "Rousettus madagascariensis"] + 2
pies_pos$y2[pies_pos$bat_species== "Rousettus madagascariensis"] <- pies_pos$latitude_s[pies_pos$bat_species== "Rousettus madagascariensis"] - 2

#pies_pos$x2[pies_pos$bat_species== "Eidolon dupreanum/Rousettus madagascariensis"] <- pies_pos$longitude_e[pies_pos$bat_species== "Eidolon dupreanum/Rousettus madagascariensis"] +3
#pies_pos$y2[pies_pos$bat_species== "Eidolon dupreanum/Rousettus madagascariensis"] <- pies_pos$latitude_s[pies_pos$bat_species== "Eidolon dupreanum/Rousettus madagascariensis"] -1

coordinate$x2[coordinate$bat_species== "Eidolon dupreanum/Rousettus madagascariensis"] <- coordinate$longitude_e[coordinate$bat_species== "Eidolon dupreanum/Rousettus madagascariensis"] +1
coordinate$y2[coordinate$bat_species== "Eidolon dupreanum/Rousettus madagascariensis"] <- coordinate$latitude_s[coordinate$bat_species== "Eidolon dupreanum/Rousettus madagascariensis"] +0.2

head(pies_pos)

#This is figure 1

#all positives
p4_pos <- p2b+
  annotate("segment", x=pies_pos$longitude_e, xend=pies_pos$x2,y=pies_pos$latitude_s,yend=pies_pos$y2,size=1)+ # put the lines
  annotate("segment", x=coordinate$longitude_e, xend=coordinate$x2,y=coordinate$latitude_s,yend=coordinate$y2,size=1)+
  geom_scatterpie(aes(x=x2, y=y2, r=(N/1000)*8), 
    data = pies_pos, cols="plot_class_pos", long_format=TRUE) +
  theme_bw() +theme(panel.grid = element_blank(),
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(.1,.5,.1,.5),"cm"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.85,.15),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.text = element_text(size = 9)) +
  scale_fill_manual(values=colz)
  # geom_scatterpie_legend(log10(c(10,200)/1.2),
  #                        x=55, y=-23.5, 
  #                        n=2,
  #                        labeller = function(x) paste(10^(x)*1,"indiv"))

p4_pos

ggsave(file = paste0(homewd, "/final_figures/supplemental/Sfig1_sampling_map_only.png"),
       plot = p4_pos,
       units="mm",
       width=60,
       height=60,
       scale=3,
       dpi=300)
