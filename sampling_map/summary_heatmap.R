rm(list=ls())

#packages
library(plyr) 
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(hrbrthemes)
library(reshape)
library(paletteer)
library(ggrepel)
library(hrbrthemes)
library(reshape)
library(paletteer)
library(ggh4x)
#####################################################################

homewd = "/Users/gwenddolenkettenburg/Desktop/developer/mada-bat-picornavirus" 

dat <- read.csv(file = paste0(homewd,"/metadata/summary_contig_reads_for_heatmap.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)

dat$genus<-factor(dat$genus, levels=c("Cardiovirus", "Hepatovirus", "Kobuvirus", "Kunsagivirus",
                                      "Mischivirus", "Sapelovirus","Sapovirus","Teschovirus","Unclassified picornavirus"))
dat$family<-factor(dat$family, levels=c("Caliciviridae", "Picornaviridae"))
dat$virus<-factor(dat$virus, levels=c("Pteropus rufus mischivirus", "Eidolon dupreanum cardiovirus", "Eidolon dupreanum hepatovirus",
                                      "Eidolon dupreanum kobuvirus", "Eidolon dupreanum kobuvirus 2", "Eidolon dupreanum kunsagivirus",
                                      "Eidolon dupreanum sapelovirus 1", "Eidolon dupreanum sapelovirus 2", "Eidolon dupreanum sapovirus 1",
                                      "Eidolon dupreanum sapovirus 2", "Eidolon dupreanum sapovirus 3", "Eidolon dupreanum sapovirus 4", "Eidolon dupreanum teschovirus 1",
                                      "Rousettus madagascariensis picornavirus 1", "Rousettus madagascariensis picornavirus 2", "Rousettus madagascariensis picornavirus 3",
                                      "Rousettus madagascariensis sapelovirus 1", "Rousettus madagascariensis sapovirus 1", "Rousettus madagascariensis sapovirus 2",
                                      "Rousettus madagascariensis sapovirus 3", "Rousettus madagascariensis teschovirus 1", "Rousettus madagascariensis teschovirus 2"))
dat$species<-factor(dat$species, levels=c("Pteropus rufus","Eidolon dupreanum", "Rousettus madagascariensis"))

ED<-subset(dat, species=="Eidolon dupreanum")
PR<-subset(dat, species=="Pteropus rufus")
RM<-subset(dat, species=="Rousettus madagascariensis")


partial<-subset(dat,type=="partial")
full<-subset(dat,type=="full")


#facet by viral family, genus, and virus to look at the data


##CONTIGS FIRST
p1_contig <- ggplot(dat, aes(x=species, y=virus, fill=num_contigs)) +
  geom_tile() +
  #geom_tile(aes(fill = cut(num_genome,breaks=0:6, labels=1:6))) +
  #scale_fill_manual(values=c("lightblue1","skyblue", "royalblue1")) +
  #scale_fill_viridis_c(option = "B", direction = -1) +
  scale_fill_gradient(low="yellow", high="red")+
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  facet_nested(family+genus~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Contigs",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="right", limits=rev)+
  scale_x_discrete(position="top")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(face="italic"),
        legend.position = "bottom")
p1_contig


#facet by viral family and genus
p2_contig <- ggplot(dat, aes(x=species, y=genus, fill=num_contigs)) +
  geom_tile() +
  #geom_tile(aes(fill = cut(num_genome,breaks=0:6, labels=1:6))) +
  #scale_fill_manual(values=c("lightblue1","skyblue", "royalblue1")) +
  #scale_fill_viridis_c(option = "B", direction = -1) +
  scale_fill_gradient(low="yellow", high="red")+
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Contigs",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="left", limits=rev)+
  scale_x_discrete(position="top")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=12,face="italic"),
        axis.text.x = element_text(size=12,face="italic"),
        legend.position = "bottom")
p2_contig


#By RPM
p1_rpm <- ggplot(dat, aes(x=species, y=virus, fill=num_reads_log)) +
  geom_tile() +
  #geom_tile(aes(fill = cut(num_genome,breaks=0:6, labels=1:6))) +
  #scale_fill_manual(values=c("lightblue1","skyblue", "royalblue1")) +
  #scale_fill_viridis_c(option = "B", direction = -1) +
  scale_fill_gradient(low="yellow", high="red")+
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  facet_nested(family+genus~., scales="free", space="free",
               switch="y", nest_line = element_line(color="white"), solo_line = TRUE)+
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Reads (log10)",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="right", limits=rev)+
  scale_x_discrete(position="top")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=7,face="italic"),
        axis.text.x = element_text(face="italic"),
        legend.position = "bottom")
p1_rpm


#facet by viral family and genus
p2_rpm <- ggplot(dat, aes(x=species, y=genus, fill=num_reads_log)) +
  geom_tile() +
  #geom_tile(aes(fill = cut(num_genome,breaks=0:6, labels=1:6))) +
  #scale_fill_manual(values=c("lightblue1","skyblue", "royalblue1")) +
  #scale_fill_viridis_c(option = "B", direction = -1) +
  scale_fill_gradient(low="yellow", high="red")+
  #facet_wrap(~family, ncol=1,  scales = "free_y")+
  # facet_nested(family~., scales="free", space="free",
  #              switch="y")+
  
  labs(#title = expression("Diversity of" ~italic(Picornavirales) ~"sequences"),
    x = "",
    y= "",
    fill="Reads (log10)",
    title="")+
  theme_linedraw()+
  scale_y_discrete(position="left", limits=rev)+
  scale_x_discrete(position="top")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12, face="italic"), 
        axis.text.y = element_text(size=12,face="italic"),
        axis.text.x = element_text(size=12,face="italic"),
        legend.position = "bottom")
p2_rpm


##Let's just show reads for now

